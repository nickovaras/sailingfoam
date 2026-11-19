# sailingfoam
OpenFOAM, local LLMs and Sailboats

Computational Fluid Dynamics (CFD) has traditionally been a field that required substantial 
effort and knowledge.

This repo hosts some experiments done to see how local LLMs can reduce the friction needed to do 
some CFD. To that purpose, I set to explore a statement I heard one of my sailing buddies utter:

_"At higher wind speeds, the jib and main start being less like a biplane (with a slot) and begin to act 
as a single airfoil"_

The experiment involved:
- A Jetson Orin (arm64)
- Qwen3.8-27B
- OpenCode
- OpenFOAM docker containers
- A Qwen-generated Python script to generate meshes

AI was used to set up the enviroment and experiments only. The CFD simulation can be run without it.
The Python script takes parameters such a angle of heel and sail trim to generate the desired mesh.


## How to Use:

On host machine:
```
docker compose up -d
docker exec -it openfoam bash
```

Within container:

```
cd cases/sailBoat
./Allclean
./Allrun
paraview
```