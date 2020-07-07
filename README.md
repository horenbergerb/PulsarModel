The goal here was/is to apply electrodynamics and visualization techniques to see what the radiation from a pulsar looks like.

The simplest model of a pulsar I could think of is derived as follows:
   A pulsar is a sphere which rotates on some axis.
   The plasma within the pulsar behaves as a current.
   The magnetic axis of the pulsar's current is offset from the rotation axis.
   Thus, I am modeling the pulsar as a rotating loop of current whose normal vector is offset from the rotation axis.

There's a lot of code here. It shows a developmental process in steps. I'd order them as follows:
  MagnetostaticStationaryLoopPotential
  ElectrodynamicyStationaryLoopPotential
  ElectrodynamicMovingLoopPotential
  CurrentLoopAnim (I actually made this first just to get a feel for the animation)
  ElectrodynamicMovingLoopPotentialAnimated
  ElectrodynamicPoyntingVectorAnimated (WIP)

TO DO:
Calculate Poynting vector and average pointing vector over a cycle
Optimize calculations (I'm keeping my eyes peeled for any mathematical or CS optimizations)
