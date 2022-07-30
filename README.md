# ANALYSIS-OF-RIGID-BODY-KINEMATICS


ANALYSIS OF RIGID BODY KINEMATICS
Simran Karim

OBJECTIVE
Kinematics is a subfield of physics that involves the motion of bodies without considering the forces that created the motions. The objective of this laboratory was to learn how to capture the relative motions between two rigid bodies and to track and analyze those relative motions using a surrogate of a mechanical knee joint.  The two rigid bodies consisted of a femur and tibia as shown in Figure 1. There are a number of ways to track the kinematics of rigid bodies. This lab utilized a 3D electromagnetic digitizer to manually record the 3D coordinates of the two rigid bodies at different positions. This, like other similar methods, is useful in the field of Biomedical Engineering as it allows identification of important information concerning spatial positions of body segments, or the rate at which the segments are moving. For example, gait analysis uses a similar tracking technique that allows physicians or researchers to understand how a disease may influence a patient’s gait. Therefore, it is an important method that aids in understanding the kinematics of rigid bodies, or in the case of this laboratory, the aforementioned technique was utilized to record measurements that were later used to find the positions and the relative angles of the two rigid bodies in a knee joint. 

METHODS
The goal of the experiment mentioned below was to establish an embedded coordinate system on the femur and the tibia with respect to the kinematic marker clusters placed on both of the rigid bodies. The experiment was repeated at three different passive flexion angles by moving the tibia relative to the fixed femur. 

	Experimental Procedure: First, the tibial segment of the mechanical knee joint was set at approximately 5˚ flexion. The flexion angle and the distances among 21 landmarks on the femur and tibia, shown in Figure 2a and 2b, were measured using a goniometer and a ruler, respectively. The landmarks consisted of the femur and tibia marker clusters shown in Figure 1, the knee joint mounting frame positions and the medial, lateral, anterior, posterior of both femur and tibia distal and proximal, depending on what was necessary. Finally, with the aid of the 3D MicroScribe digitizer, those landmarks were manually digitized in a specific order (Figure 2a and 2b) by directly placing the MircroScribe probe on the surface of the femur and tibia, and the data of the spatial positions of those body segments were saved in a computer for later analysis. The exact procedure was repeated for flexion angles at ~30˚ and ~90˚ by only lowering the tibial segment of the knee joint.













	Data Analysis Overview: The goal of this data analysis was to first calculate the transformation matrix between the femoral and tibial marker clusters Tfm→tm(t) as well as the transformation matrix between the embedded coordinate systems between the two rigid bodies Tfe→fm(t) and Ttm→te(t), shown in Figure 4. In order to do that, the programming language MATLAB was utilized. First, all the digitized points for position 1, 2 and 3 were loaded in a MATLAB script named Euler_Digitizer_2021_Karim_Simran.m. Using those points, two vectors in the direction of the axis shown in Figure 2a, 2b and 3 were calculated. Taking cross product of those two vectors yielded in one vector that was perpendicular to the previous vectors. The new vector was used to calculate the 3rd perpendicular vector in the 2nd cross product. The resulting vectors in the x,y,z axis were then normalized to create an orthonormal coordinate system for the femur marker set, tibia marker set, mounting frame corners, embedded distal femur and embedded proximal tibia. Using those normalized data points, the three transformation matrices Tfm→tm(t), Tfe→fm(t) and Ttm→te(t) were calculated with the function transform_matrix.m. Once the three important transformation matrices were obtained, the final part of the analysis included finding the transformation matrix between the femoral and tibial embedded systems Tfe→te(t) and the Euler angles associated with it. Equation 1 was utilized to calculate the transformation matrix Tfe→te(t). Finally, using the Tfe→te(t) transformation matrix the Euler angles were calculated for a counter-clockwise rotation, using the convention shown in Equation 2. The Euler angles results were plotted in MATLAB for the three different positions. Additionally, the magnitude of each translation vector in the Tfe→te(t) transformation matrix were also obtained using the 4th column of the matrix. 

                    Equation 1

							        Equation 2
RESULTS

 
Plot 1.  Representation of the Euler angles Ψ, θ, ϕ in plot a, b and c, respectively, for position 1,2,3.

Table 1. Translation Vector Magnitudes   
	Magnitude of Translation Vector [𝐓𝐟𝐞→𝐭𝐞] (mm) 
Group 4 Data Position 1 	153.4955
Group 4 Data Position 2 	143.5246
Group 4 Data Position 3 	119.6461


Table 2. Manual Calculations of Euler Angles

	Position 1	Position 2	Position 3
Ψ angles [degrees]	10.1379	24.6135	93.0340
Φ angles [degrees]	0.4435	18.2148	45.0229




DISCUSSSION
In order to compare the MATLAB analysis of the Euler angles, the Euler angles for 2 different rotations were manually calculated, as shown in Table 2. First, the subplot Plot a in the figure Plot 1, obtained using MATLAB, shows values of the Ψ angle ~ -9.2544˚, ~ -23.639˚and ~ -90.5794˚ for position 1, 2 and 3, respectively. Comparing this data with the manually calculated Ψ angles for all positions shown in Table 2, it was calculated that the Ψ angle values are slightly off by a factor of ~0.95, which shows that the two groups of data are close to each other. Comparing these two groups of data with the measured data obtained using a goniometer for the flexion angle between the femur and tibia (0, 30˚, 90˚ for position 1, 2, 3 respectively), the values are off more than a factor ~0.95 but not significantly for the data to be incorrect. Using these interpretations, it can be said that the measurement of the flexion angles using the goniometer is not an accurate depiction of the flexion angles 0˚,30˚,90˚due to human error. The MATLAB data and the manually obtained data are more reliable sources to understand thee spatial positions and angles of the femur and tibial rotations since the data between these two groups are close to each other. The MATLAB values of Φ angle shown in subplot Plot b of Plot 1 were ~ 1.5036˚, ~ 15.8380˚ and ~ 42.3199 for position 1, 2 and 3, respectively. Comparing this data with the manually calculated Φ angles for all positions shown in Table 2, it was calculated that the Φ angle values are slightly off by a factor of ~0.7, which can be used to interpret that the two groups of data are close to each other. Finally, comparing these two groups of data with the measured data for tibial rotation with was obtained with the aid of a goniometer (~0˚, ~3˚, ~45˚ for position 1, 2, 3 respectively), it is noticeable that the measured data and the two groups of calculated data are off significantly for position 2, but not so much for position 1 and 3. Using these interpretations, it can be conveyed that goniometer method to obtain the angles may have potentially given inaccurate values since the MATLAB and the manually calculated data were close to each other, therefore, making the latter two methods  reliable.

The translation vector magnitude calculated using MATLAB shown in Table 1 gave an idea of the femur and tibial embedded coordinate system’s origin.

The mounting frame data was also obtained, even though not analyzed, since the mounting frame data provided a constant coordinate system for the mechanical knee movements. In order words, it served a control for the experiment. This allowed for a reliable experiment to obtain close-to-accurate data. 

In conclusion, the data obtained using MATLAB seem very reliable. However, there is one potential for unreliability on these MATLAB calculated data since true values were not accurately obtained using the goniometer. As a result, it is unsure to tell if the MATLAB calculated reflect the true values regarding the position and the angles at the different positions. The same can be said for the translated vector magnitude data. Even though the magnitude of the translation vector was obtained, it was not possible to compare these with the true value of the distances since there are potential sources of human error when using a ruler.  One potential solution to this issue would be to use a digital goniometer/ ruler or a similar device that allows physically measuring the rotational angles and distances between points with high accuracy. That way future experimenters will potentially have the true values to use as comparisons for their MATLAB calculated data. 





![image](https://user-images.githubusercontent.com/105514187/181907641-23ba9d68-e8a1-4bc1-84a6-af0154b59218.png)
