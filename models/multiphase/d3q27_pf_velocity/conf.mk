ADJOINT=0
TEST=FALSE
OPT="(OutFlow+BGK+thermo*planarBenchmark+ferro)*autosym"
# SC: Solid Contact
# 	This option currently fixes the bottom layer of nodes to be 
# 	solid with the contact angle defined in input.
# thermo: thermocapillary flows
# 	Options resolves the temperature field with an RK4 integration
# 	and updates the surface tension as a result
# ferro: ferrofluid flows
# 	Option to resolve magnetic field and its impact on magnetisable
# 	droplets
