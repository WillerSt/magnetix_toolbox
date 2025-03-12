import fmt_scen as scen
import matplotlib.pyplot as plt
import sys

results = scen.result_xml(sys.argv[1])
startIdx = 100
Bquery = results.gather_point_queries('B')
Bcomp = scen.split_components(Bquery[0])
BxP34 = Bcomp[0]
ByP34 = Bcomp[1]

Bcomp = scen.split_components(Bquery[3])
BxP12 = Bcomp[0]
ByP12 = Bcomp[1]

measurements = scen.result_xml('input/validation/TEAM_Problem_32_case_3_flux.xml')

Bx_c1c2 = measurements.get_scalar_quantity('Bx_c1c2')
By_c1c2 = measurements.get_scalar_quantity('By_c1c2')
Bx_c3c4 = measurements.get_scalar_quantity('Bx_c3c4')
By_c3c4 = measurements.get_scalar_quantity('By_c3c4')


plt.figure()
plt.plot()
plt.plot(BxP34[startIdx:],ByP34[startIdx:])
plt.plot(Bx_c3c4, By_c3c4,color='red')
plt.show(block=False)
plt.axis('equal')
plt.xlabel('$B_x$ in T')
plt.ylabel('$B_y$ in T')
plt.show(block=False)

plt.figure()
plt.plot()
plt.plot(BxP12[startIdx:],ByP12[startIdx:])
plt.plot(Bx_c1c2, By_c1c2,color='red')
plt.axis('equal')
plt.xlabel('$B_x$ in T')
plt.ylabel('$B_y$ in T')
plt.show()



