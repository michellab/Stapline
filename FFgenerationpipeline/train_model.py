import numpy as np
import openmm as mm
import parmed
import torch
import torch.nn as nn
import torch.optim as optim


class PositiveFloatConstraint:
    def __call__(self, parameter):
        return torch.clamp(parameter, min=0.0)


class IntegerConstraint:
    def __call__(self, parameter):
        return torch.round(parameter)


class PhaseConstraint:
    def __call__(self, parameter):
        return 180.0 * torch.round(parameter / 180.0)


class StructureModifier(torch.nn.Module):
    def __init__(self, structure: parmed.Structure):
        super().__init__()
        self.structure = structure
        self.phi_k = nn.Parameter(
            torch.tensor([dih.phi_k for dih in self.structure.dihedral_types]),
            constraint=PositiveFloatConstraint(),
        )
        self.per = nn.Parameter(
            torch.tensor([dih.per for dih in self.structure.dihedral_types]),
            constraint=IntegerConstraint(),
        )
        self.phase = nn.Parameter(
            torch.tensor([dih.phase for dih in self.structure.dihedral_types]),
            constraint=PhaseConstraint(),
        )

    def forward(
        self,
        xyz,
    ):
        self.update_structure()
        current_MM_energies = self.compute_energy(xyz)
        return current_MM_energies
        # Update dihedral parameters

    def update_structure(self):
        for dih, phi_k, per, phase in self.structure.dihedral_types:
            dih.phi_k = phi_k.item()
            dih.per = per.item()
            dih.phase = phase.item()

    def compute_MM_energy(self, frames):
        parm = self.structure
        system = parm.createSystem(nonbondedMethod=None)

        # Make the context and set the positions
        context = mm.Context(system, mm.VerletIntegrator(0.001))

        energy_total = []
        for positions in frames:
            context.setPositions(positions)
            # Find the energy decomposition
            energy_total.append(
                parmed.openmm.energy_decomposition(parm, context)["total"]
            )
        return energy_total

    def writefrcmod(self, file_out="test.frcmod"):
        parmed.tools.writeFrcmod(self.structure, file_out)


structure_modifier = StructureModifier(structure)
optimizer = torch.optim.Adam(structure_modifier.parameters(), lr=learning_rate)

for epoch in range(num_epochs):
    for i, (positions, target_MM_energies) in enumerate(train_loader):
        optimizer.zero_grad()
        predicted_MM_energies = structure_modifier.compute_MM_energy(positions)
        loss = F.mse_loss(predicted_MM_energies, target_MM_energies)
        loss.backward()
        optimizer.step()


import torch.optim as optim

# Define your input and target data
input_data = ...  # list of energy values that you want to fit to
target_data = ...  # list of corresponding atom positions

# Define your model and optimizer
model = StructureModifier(structure)
optimizer = optim.Adam(model.parameters(), lr=0.001)

# Define your loss function
loss_fn = nn.MSELoss()

# Train the model
for epoch in range(num_epochs):
    running_loss = 0.0
    for i in range(len(input_data)):
        # Zero the parameter gradients
        optimizer.zero_grad()

        # Forward pass
        output = model(target_data[i])

        # Compute the loss
        loss = loss_fn(output, input_data[i])

        # Backward pass and optimize
        loss.backward()
        optimizer.step()

        # Print statistics
        running_loss += loss.item()
        if i % 100 == 99:
            print(f"Epoch {epoch+1}, batch {i+1}: loss={running_loss/100:.4f}")
            running_loss = 0.0


structure_modifier.writefrcmod(file_out="optimised dihedrals.frcmod")
