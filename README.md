# Determining the limits of stability of the circular restricted three-body problem
This is my research project that I did at the Loudoun Academy of Science. The code is all written in IDL and requires the IDL IDE to run.

This is a two dimensional simulation of a binary system with a relatively massless satellite.

<h2>Abstract</h2>
The elliptical restricted three-body problem describes a system of two heavy objects, with a lighter one in orbit around one of the heavy objects. The two heavy objects have coplanar, elliptical orbits and are unaffected by the gravitational field of the lighter object. The problem at hand was to determine the relationship between system stability and the mass ratio of the two heavy objects for systems in which the smallest object orbits around one of the other two masses. This was pursued through the programming of both retrograde and prograde simulations using Interactive Data Language (IDL). A fourth order Runge-Kutta method was implemented to estimate the motions of the objects. An adaptive-stepsize control routine was also added to optimize computational time and accuracy. The system was tested multiple times with different initial conditions for the orbital radius of the lightest object. It was discovered that retrograde orbits are always more stable than their corresponding prograde orbits, and systems in which the satellite orbits the smaller mass are less stable than those where the satellite orbits the larger mass. These results can assist in determining the realistic limits of the three-body problem, which can provide guidelines for astrophysicists in creating realistic simulations relevant to binary systems with satellites.<br>
<h2>Example Image 1</h2>
Binary ass ratio: 100:1<br/>
Orbital radius of <i>m3</i> over orbital radius of <i>M1</i>: 0.1<br/>
Orbital direction: retrograde<br/>
Coordinate system: synodic<br/>
<img src="m100v135r1rs copy.png"/><br/>
<h2>Example Image 2</h2>
Binary ass ratio: 1:1<br/>
Orbital radius of <i>m3</i> over orbital radius of <i>M1</i>: 0.7<br/>
Orbital direction: retrograde<br/>
Coordinate system: synodic<br/>
<img src="m1v15r7rs copy.png"/>
