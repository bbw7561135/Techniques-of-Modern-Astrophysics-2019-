
import matplotlib.pyplot as plt
import numpy as np

#-----------------------------------------------

# Files

filelist=[]
for i in range(0,20):
    filelist.append("star%s.txt" %i)


mesafiles = [r'/home/samantha/4910U/Project1/Profiles/hist1.data',r'/home/samantha/4910U/Project1/Profiles/hist2.data',r'/home/samantha/4910U/Project1/Profiles/hist3.data',r'/home/samantha/4910U/Project1/Profiles/hist4.data',r'/home/samantha/4910U/Project1/Profiles/hist6.data', r'/home/samantha/4910U/Project1/Profiles/hist7.data',r'/home/samantha/4910U/Project1/Profiles/hist8.data',r'/home/samantha/4910U/Project1/Profiles/hist9.data',r'/home/samantha/4910U/Project1/Profiles/hist10.data',r'/home/samantha/4910U/Project1/Profiles/hist11.data',r'/home/samantha/4910U/Project1/Profiles/hist12.data',r'/home/samantha/4910U/Project1/Profiles/hist13.data',r'/home/samantha/4910U/Project1/Profiles/hist14.data',r'/home/samantha/4910U/Project1/Profiles/hist15.data',r'/home/samantha/4910U/Project1/Profiles/hist16.data',r'/home/samantha/4910U/Project1/Profiles/hist17.data',r'/home/samantha/4910U/Project1/Profiles/hist19.data',r'/home/samantha/4910U/Project1/Profiles/hist20.data',r'/home/samantha/4910U/Project1/Profiles/hist21.data',r'/home/samantha/4910U/Project1/Profiles/hist22.data'] 


num = len(mesafiles)

models = [783,891,792,904,942,807,911,670,586,772,765,857,982,776,762,947,925,980,914,903]
 
for i in range(0,num): 
	
	data = np.genfromtxt(mesafiles[i], skip_header=6)  
]
	star_mass = data[models[i],2]
	log_L = data[models[i],3]
	log_R = data[models[i],4]
	log_Teff = data[models[i],5]
	log_Tc = data[models[i],6] 
	log_P = data[models[i],8]
	X_H = data[models[i],9]
	X_He3 = data[models[i],10]
	X_He4 = data[models[i],11]


	mesa=plt.plot(log_Teff,log_L,'bo',label='MESA'))
	plt.xlabel('log $T_{eff}$ [kK]')
	plt.ylabel('log(L/$L_{\odot}$)')
	plt.title('H-R Diagram Comparing MESA and Our Models')

	

for fname in filelist:
	data=np.loadtxt(fname, unpack = True)
	X=data[0]
	Y=data[1]
	ourmodels=plt.plot(np.log10(X),np.log10(Y),'ro',label='Our Stars')


plt.legend([mesa[0],ourmodels[0]],['MESA', 'Our Stars'])
plt.gca().invert_xaxis()

plt.savefig(r'/home/samantha/4910U/Project1/proj1_hr_both.pdf')

plt.show() 




