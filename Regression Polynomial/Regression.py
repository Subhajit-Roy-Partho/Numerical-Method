import math
import os
x=[1.2,2.1,2.8,4.1,4.9,6.2,7.1,7.9,8.9]
y=[4.2,6.8,9.8,13.4,15.5,19.6,21.6,25.4,28.6]


# print(len(x),len(y))
def sigma(x,n):
	p=0
	for m in range(0,len(y)):
		p+=(x[m])**(n)
	return p
# print(sigma(y,1))

def sigma2(x,n,y,m):
	p=0
	for j in range(0,len(y)):
		p+=math.pow(x[j],n)*math.pow(y[j],m)
		# p+=(((x[j])**(n)))*(((y[j])**[m]))
	return p
# print(sigma2(x,0,y,1))
print("Enter the number of degrees")
# deg=input()
deg=2

a=[]
b=[]
for i in range(0,deg):
	b.append(sigma(x,i))
b.append(sigma2(y,1,x,0))
# print(b)
a.append(b)




for pp in range(0,deg):
	b=[]
	for i in range(1+pp,deg+1+pp):
		b.append(sigma(x,i))
	b.append(sigma2(y,1,x,pp+1))
	# print(b)
	a.append(b)
# print(a[:])




def interows(array,x,y):
	b=array[y][:]
	array[y][:]=array[x][:]
	array[x][:]=b

array=a
print("Enter number of degrees")
# deg=input()
# deg=3
equ=deg

# file=open("input.txt","r")
# for i in range(0,deg):
# 	file.seek(0,0)
# 	line=file.readlines()[i]
# 	store=line.split()
# 	length=len(store)
# 	for k in range(0,length):
# 		store[k]=float(store[k])
# 	array.append(store)


for m in range(1,deg):


	maxi=[]
	for i in range(m-1,deg):
		maxi.append(array[i][m-1])
	huge=maxi.index(max(maxi))+(m-1)
	interows(array,m-1,huge)
	del maxi


	for p in range(m,deg):
		const=(array[p][0+m-1] / array[0+m-1][0+m-1])
		# print(const)
		
		for o in range(0,equ+1):
			array[p][o]-=(const*array[0+m-1][o])
			
		# print(array)

result=[]
count=0
for n in range(deg-1,-1,-1):
	b=array[n][deg]
	# print(b)
	count=0
	for p in range(deg-1,n,-1):		
		b-=array[n][p]*result[count]
		count+=1
	net=b/array[n][n]
	result.append(net)
	

result.reverse()
print(str("X=\n")+str(result))

file=open("plot2.gnplt","w")
file.seek(0,0)
file.write(str('set terminal postscript\nset output "Out.ps"\nset xlabel "X"\nset ylabel "Y"\nplot "input.dat" using 1:2\nreplot '))
for h in range(0,deg):
	if h == 0:
		file.write(str(result[h])+str(' + '))
	if h == deg-1:
		file.write(str(result[h])+str('*x**')+str(h))
	if h != 0 and h !=deg-1:
		file.write(str(result[h])+str('*x**')+str(h)+str(' + '))

file.write('\nquit')
file.close()
os.system("gnuplot plot2.gnplt")
os.system("ps2pdf Out.ps out.pdf")
os.system("rm plot2.gnplt Out.ps")
print("Look out for 'out.pdf'")

#with degree 5 intersecting almost all the points
#X=
#[-0.7097936010778282, 4.320529344814003, -0.2987872619743378, 0.01858070804831512, 0.0002466319480689103]
#with degree 2 
#Intersecting only the first point.
# X=
# [0.5063551688451235, 3.104929280539688]

# Slop of the graph = 3.104929280539688
#

#Gnuplot values
#slop =3.10493
#x = 0.506355, 3.10493