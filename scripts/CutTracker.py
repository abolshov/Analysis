import uproot
import numpy as np
import vector
import awkward as ak
import matplotlib.pyplot as plt
import os


sample_type_map = {1: "GluGluToRadion",
                   2: "GluGluToBulkGraviton", 
                   3 : "VBFToRadion", 
                   4 : "VBFToBulkGraviton"}


class CutTracker():
    def __init__(self, file_name):
        file = uproot.open(file_name)
        self.tree = file["Events"]
        self.branches = self.tree.arrays()
        self.cutflow = {"total": len(self.branches)}
        self.masspoint = int(self.branches['X_mass'][0]) if 'X_mass' in self.tree.keys() else None
        self.sample_key = self.branches['sample_type'][0]
        self.sample_type = sample_type_map[self.sample_key] if self.sample_key in sample_type_map else None
        
    def Apply(self, func, cut_name):
        if cut_name not in self.cutflow:
            mask = func(self.branches)
            self.branches = self.branches[mask]
            self.cutflow[cut_name] = len(self.branches)
        else:
            print(f"Cut {cut_name} already in cutflow table, no action taken")
            
    def Print(self):
        cut_lst = list(self.cutflow.keys())
        last_cut = cut_lst[-1]
        print(f"{last_cut}: {self.cutflow[last_cut]}")

    def GetLastCutDescription(self):
        cut_lst = list(self.cutflow.keys())
        return cut_lst[-1]

    def PlotCutflow(self, show=False, percentages=False, save_path=None):
        cut_names = self.cutflow.keys()
        passed_evt_info = list(self.cutflow.values())
        if percentages:
            passed_evt_info = [val/passed_evt_info[0]*100 for val in passed_evt_info]
        fig, ax = plt.subplots()
        bar_container = ax.bar(cut_names, passed_evt_info)
        title = f"Cutflow {self.sample_type} M={self.masspoint}" if (self.sample_type and self.masspoint) else "Cutflow"
        ax.set(ylabel='Number of events passed', title=title)
        if percentages:
            ax.bar_label(bar_container, fmt=lambda x: f"{x:.1f}%", fontsize=8)
        else:
            ax.bar_label(bar_container, fmt=lambda x: int(x), fontsize=8)
        plt.xticks(rotation=45, ha='right')

        if show:
            plt.show()

        if self.sample_type and self.masspoint:
            if save_path:
                plt.savefig(os.path.join(save_path, f"cutflow_{self.sample_type}_M{self.masspoint}.pdf"), format='pdf', bbox_inches='tight')
            else:
                plt.savefig(f"cutflow_{self.sample_type}_M{self.masspoint}.pdf", format='pdf', bbox_inches='tight')
        else:
            if save_path:
                plt.savefig(os.path.join(save_path, f"{save_path}cutflow_{self.sample_key}.pdf"), format='pdf', bbox_inches='tight')
            else:
                plt.savefig(f"cutflow_{self.sample_key}.pdf", format='pdf', bbox_inches='tight')
        plt.close()

    def GetSelectedEvents(self):
        return self.branches

    def NumEventsSelected(self):
        return len(self.branches)