# Nuetron-Path-Simulation

Made use of random number generators in Python to simulate the three-dimensional random path of a thermal neutron in the wall of a nuclear reactor. This could be used to test the effectiveness of different materials when absorbing a thermal neutron.
Made use of Monte Carlo integration to create an exponential distribution to choose the length at which a neutron would move between each collision.
Inserted the quantum cross-section of each material into the exponential distribution to calculate the proportion of neutrons transmitted, absorbed or scattered. Calculated the characteristic attenuation length with error analysis for different materials.
