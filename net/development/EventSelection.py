import awkward as ak
import numpy as np
from typing import Dict, Any
from Utils import dr2

class Run3bbWWSelector:
    """
    Class for performing analysis event selections for Run3 X->HH->bbWW analysis.
    """
    def __init__(self,
                 channel: str,
                 cfg: Dict[str, Any]):
        self.channel = channel
        self.cfg = cfg

    def _apply_gen_level(self,
                         chunk: ak.Array) -> ak.Array:
        gen_cfg = self.cfg['gen_level']
        
        mask = (chunk['genb1_pt'] > gen_cfg['bquark_pt']) & \
               (np.abs(chunk['genb1_eta']) < gen_cfg['bquark_eta']) & \
               (chunk['genb2_pt'] > gen_cfg['bquark_pt']) & \
               (np.abs(chunk['genb2_eta']) < gen_cfg['bquark_eta'])

        thresh2 = gen_cfg['boosted_thresh']**2
        match gen_cfg['bb_topology']:
            case 'boosted':
                mask &= (dr2(data=chunk, prefix1='genb1', prefix2='genb2') < thresh2)
            case 'resolved':
                mask &= (dr2(data=chunk, prefix1='genb1', prefix2='genb2') >= thresh2)
            case 'mixed':
                pass
            case _:
                raise ValueError(f'Illegal `bb_topology` value: {gen_cfg['bb_topology']}.')

        match self.channel:
            case 'DL':
                # only left to check leptons here
                mask &= (
                    (chunk['genV1prod1_pt'] > gen_cfg['lep_pt']) & (chunk['genV2prod1_pt'] > gen_cfg['lep_pt'])
                )
            case 'SL':
                # here need to check lepton first 
                mask &= (chunk['genV1prod1_pt'] > gen_cfg['lep_pt'])

                # and now check W->qq quarks
                match gen_cfg['qq_topology']:
                    case 'boosted':
                        mask &= (dr2(data=chunk, prefix1='genV2prod1', prefix2='genV2prod2') < thresh2)
                    case 'resolved':
                        mask &= (dr2(data=chunk, prefix1='genV2prod1', prefix2='genV2prod2') >= thresh2)
                    case 'mixed':
                        pass
                    case _:
                        raise ValueError(f'Illegal `qq_topology` value: {gen_cfg['qq_topology']}')
                    
        return chunk[mask]

    def _apply_reco_level(self,
                          chunk: ak.Array) -> ak.Array:
        reco_cfg = self.cfg['reco_level']

        num_ak4_bjets = ak.sum(
            (chunk['centralJet_btagPNetB'] > reco_cfg['ak4_btag']) &
            (chunk['centralJet_pt'] > reco_cfg['ak4_pt']) &
            (np.abs(chunk['centralJet_eta']) > reco_cfg['ak4_b_eta']),
            axis=1
        )

        num_ak8_bjets = ak.sum(
            (chunk['SelectedFatJet_btagPNetB'] > reco_cfg['ak8_btag']) &
            (chunk['SelectedFatJet_pt'] > reco_cfg['ak8_pt']) &
            (np.abs(chunk['SelectedFatJet_eta']) > reco_cfg['ak8_b_eta']),
            axis=1
        )

        if self.channel == 'DL':
            mask = (num_ak4_bjets >= 2) | (num_ak8_bjets >= 1)
            mask &= (
                (chunk['lep1_pt'] > reco_cfg['lep_pt']) & (chunk['lep2_pt'] > reco_cfg['lep_pt'])
            )
        elif self.channel == 'SL':
            num_ak4_Wjets = ak.sum(
                (chunk['centralJet_pt'] > reco_cfg['ak4_pt']) &
                (np.abs(chunk['centralJet_eta']) > reco_cfg['ak4_light_eta']),
                axis=1
            )

            num_ak8_Wjets = ak.sum(
                (chunk['SelectedFatJet_pt'] > reco_cfg['ak8_pt']) &
                (np.abs(chunk['SelectedFatJet_eta']) > reco_cfg['ak8_light_eta']),
                axis=1
            )
        
            mask = (
                ((num_ak4_bjets >= 2) & (num_ak4_Wjets >= 2)) | 
                ((num_ak8_bjets >= 1) & (num_ak4_Wjets >= 2)) |
                ((num_ak4_bjets >= 2) & (num_ak8_Wjets >= 1)) |
                ((num_ak8_bjets >= 1) & (num_ak8_Wjets >= 1))
            )

            mask &= (chunk['lep1_pt'] > reco_cfg['lep_pt'])
        else:
            raise ValueError(f'Illegal channel {self.channel}.')

        return chunk[mask]

    def _apply_gen_reco_matching(self,
                                 chunk: ak.Array) -> ak.Array:
        matching_cfg = self.cfg['gen_reco_matching']

        # lepton matching
        mask = (chunk['lep1_legType'] == chunk['lep1_gen_kind'])
        if self.channel == 'DL':
            mask &= (chunk['lep2_legType'] == chunk['lep2_gen_kind'])

        ak4_thresh2 = matching_cfg['match_ak4_thresh']**2
        ak8_thresh2 = matching_cfg['match_ak8_thresh']**2

        active_chunk = chunk[mask]
        if len(active_chunk) == 0:
            return active_chunk
        
        def get_min_info(data: ak.Array, 
                         p1: str, 
                         p2: str) -> tuple:
            """Helper to get both min dR^2 and index in one pass."""
            dists = dr2(data=data, prefix1=p1, prefix2=p2)
            # We compute this once and store it to avoid re-calculating dr2
            return ak.min(dists, axis=1), ak.argmin(dists, axis=1) 
    
        b1_min_dr2, b1_idx = get_min_info(active_chunk, 'genb1', 'centralJet')
        b2_min_dr2, b2_idx = get_min_info(active_chunk, 'genb2', 'centralJet')
        fatH_min_dr2, fatH_idx = get_min_info(active_chunk, 'genHbb', 'centralJet')

        if self.channel == 'DL':
            match_mask = (
                ((b1_idx != b2_idx) & (b1_min_dr2 < ak4_thresh2) & (b2_min_dr2 < ak4_thresh2)) | 
                (fatH_min_dr2 < ak8_thresh2)
            )
        elif self.channel == 'SL':
            q1_min_dr2, q1_idx = get_min_info(active_chunk, 'genV2prod1', 'centralJet')
            q2_min_dr2, q2_idx = get_min_info(active_chunk, 'genV2prod2', 'centralJet')
            fatW_min_dr2, fatW_idx = get_min_info(active_chunk, 'genV2', 'centralJet')

            match_mask = (
                (
                    (b1_idx != b2_idx) & (q1_idx != q2_idx) & 
                    (b1_idx != q1_idx) & (b1_idx != q2_idx) & 
                    (b2_idx != q1_idx) & (b2_idx != q2_idx) &
                    (b1_min_dr2 < ak4_thresh2) & (b2_min_dr2 < ak4_thresh2) &
                    (q1_min_dr2 < ak4_thresh2) & (q2_min_dr2 < ak4_thresh2)
                ) |
                (
                    (b1_idx != b2_idx) & (b1_min_dr2 < ak4_thresh2) & (b2_min_dr2 < ak4_thresh2) & (fatW_min_dr2 < ak8_thresh2)
                ) | 
                (
                    (q1_idx != q2_idx) & (q1_min_dr2 < ak4_thresh2) & (q2_min_dr2 < ak4_thresh2) & (fatH_min_dr2 < ak8_thresh2)
                ) |
                (
                    (fatW_idx != fatH_idx) & (fatH_min_dr2 < ak8_thresh2) & (fatW_min_dr2 < ak8_thresh2)
                )
            )
        else:
            raise ValueError(f'Illegal channel {self.channel}.')
        
        return active_chunk[match_mask]

    def apply_selections(self, 
                         chunk: ak.Array) -> ak.Array:
        if self.cfg.get('gen_level'):
            chunk = self._apply_gen_level(chunk)
        if self.cfg.get('reco_level'):
            chunk = self._apply_reco_level(chunk)
        if self.cfg.get('gen_level') and self.cfg.get('reco_level') and self.cfg.get('gen_reco_matching'):
            chunk = self._apply_gen_reco_matching(chunk)
        return chunk