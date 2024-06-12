# Build llm agents that can run MD simulations 
Given a uniprot ID, these agents should be able to run a simulation. First iteration will consider all things being equal, that is perfect protein with pdb structures etc. Other iterations will include checking for alphafold structures where available. User just has to specify uniprot ID and length of production run. 

# Initial list of remaining tasks
1. Prepare protein using OSS
2. Embed in membrane or solution box (OpenMM due to python compatibility?)
3. Run standard equilibration (allow user to change)
4. Production run

# List 2
1. Check for Alphafold structures
2. Fill in missing loops/residues if indicated
3. Analysis modules (MDtraj, MDAnalysis)
