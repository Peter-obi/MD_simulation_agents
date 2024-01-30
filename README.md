# Build llm agents
Given a uniprot ID, these agents should be able to run a simulation. First iteration will consider all things being equal, that is perfect protein with pdb structures etc. Other iterations will include checking for alphafold structures where available. User just has to specify length of uniprot ID and length of production run. 

# List of tasks, will update as I finish each one
1. Is protein soluble or membrane bound? Use appropriate extraction scripts.
2. Choose pdb with highest resolution from returned results and prepare.
