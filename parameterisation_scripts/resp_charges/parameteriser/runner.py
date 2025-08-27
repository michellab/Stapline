import os
import rdkit
import psiresp

from parameteriser.prep import build_molecule_from_smiles


class FF_Generator:
    def __init__(self, input_string):

        self.input_string: str = input_string

        mol, _, backbone_list, capping_list, sidechain_list = build_molecule_from_smiles(input_string)
        self.mol: rdkit.Molecule = mol
        self.backbone_list: list = backbone_list
        self.capping_list: list = capping_list
        self.sidechain_list: list = sidechain_list
        self.psiresp_job: psiresp.job.Job

        self.charges: list


    def get_resp_job(self): 
        psirespmol = psiresp.Molecule.from_rdkit(self.mol)
        constraints = psiresp.ChargeConstraintOptions(symmetric_atoms_are_equivalent=True)
        constraints.add_charge_sum_constraint_for_molecule(psirespmol, charge=0, indices=self.capping_list)

        psirespmol.optimize_geometry=True

        geometry_options = psiresp.QMGeometryOptimizationOptions(
        method="hf",
        basis="6-31g*")

        esp_options = psiresp.QMEnergyOptions(
        method="hf",
        basis="6-31g*")

        job_multi = psiresp.Job(
            molecules=[psirespmol],
            charge_constraints=constraints,
            qm_optimization_options=geometry_options,
            qm_esp_options=esp_options,
            working_directory="resp_qm_calculations")
        
        self.psiresp_job = job_multi

    def run_resp_job(self, resp_job):
        resp_job.generate_conformers()
        resp_job.generate_orientations()
        resp_job.run()

    def run_resp_optimisation(self, resp_job):
        os.chdir("resp_qm_calculations/optimization/")
        os.system("bash run_optimization.sh")
        os.chdir("../../")
        resp_job.run()
    
    def run_qm_get_charges(self, resp_job):
        os.chdir("resp_qm_calculations/single_point/")
        os.system("bash run_single_point.sh")
        os.chdir("../../")
        resp_job.run()
        self.charges = resp_job.molecules[0].stage_2_restrained_charges
        


