#Dependencies: None
"""
/***************************************************************************
	Author			:Charles B. Cossé

	Website			:http://netdispenser.github.io

	Email			:ccosse (at) gmail.com

	Copyright		:(C) 2010-2017 Asymptopia Software.

	License			:Apache 2.0

 ***************************************************************************/
"""

import copy
import random

N=64
visc=0.0
diff=0.0
dt=0.1

def print_array(x):
	for i in range(1,N+1):
		for j in range(1,N+1):
			#print "(%d,%d) "%(i,j),
			print("%f "%x[IX(i,j)]),
		print('')
	print('')

def IX(i,j):
	return ((i)+(N+2)*(j))

def SWAP(x0,x):
	tmp=x0
	x0=x
	x=tmp
	return x0,x

def set_bnd(b,x):
	for i in range(1,N+1):
		if b==0:
			x[IX(0  ,i)] = x[IX(1,i)]
			x[IX(N+1,i)] = x[IX(N,i)]
			x[IX(i,0  )] = x[IX(i,1)]
			x[IX(i,N+1)] = x[IX(i,N)]
		elif b==1:
			x[IX(0  ,i)]=-x[IX(1,i)]
			x[IX(N+1,i)]=-x[IX(N,i)]
			x[IX(i,0  )]=x[IX(i,1)]
			x[IX(i,N+1)]=x[IX(i,N)]
		elif b==2:
			x[IX(0  ,i)]=x[IX(1,i)]
			x[IX(N+1,i)]=x[IX(N,i)]
			x[IX(i,0  )]=-x[IX(i,1)]
			x[IX(i,N+1)]=-x[IX(i,N)]

	x[IX(0  ,0  )] = 0.5*(x[IX(1,0  )]+x[IX(0  ,1)]);
	x[IX(0  ,N+1)] = 0.5*(x[IX(1,N+1)]+x[IX(0  ,N)]);
	x[IX(N+1,0  )] = 0.5*(x[IX(N,0  )]+x[IX(N+1,1)]);
	x[IX(N+1,N+1)] = 0.5*(x[IX(N,N+1)]+x[IX(N+1,N)]);

def add_source(x,s):
	size=(N+2)*(N+2)

	for i in range(size):
		x[i] += dt*s[i];


def lin_solve(b,x,x0,a,c):
	for k in range(20):
		for i in range(1,N+1):
			for j in range(1,N+1):
				x[IX(i,j)] = (x0[IX(i,j)] + a*(x[IX(i-1,j)]+x[IX(i+1,j)]+x[IX(i,j-1)]+x[IX(i,j+1)]))/c;

		set_bnd ( b, x );

def diffuse(b,x,x0):
	a=dt*diff*N*N;
	lin_solve ( b, x, x0, a, 1+4*a )

def advect(b,d,d0,u,v):
	dt0 = dt*N
	for i in range(1,N+1):
		for j in range(1,N+1):

			x=i-dt0*u[IX(i,j)]
			y=j-dt0*v[IX(i,j)]

			if x<0.5:x=0.5
			if x>N+0.5:x=N+0.5
			i0=int(x)
			i1=i0+1

			if y<0.5:y=0.5
			if y>N+0.5:y=N+0.5
			j0=int(y)
			j1=j0+1

			s1=x-i0
			s0=1-s1
			t1=y-j0
			t0=1-t1

			d[IX(i,j)]=s0*(t0*d0[IX(i0,j0)]+t1*d0[IX(i0,j1)])+s1*(t0*d0[IX(i1,j0)]+t1*d0[IX(i1,j1)])

	set_bnd ( b, d )

def project(u,v,p,div):
	for i in range(1,N+1):
		for j in range(1,N+1):
			div[IX(i,j)] = -0.5*(u[IX(i+1,j)]-u[IX(i-1,j)]+v[IX(i,j+1)]-v[IX(i,j-1)])/N;
			p[IX(i,j)] = 0.;

	set_bnd( 0, div )
	set_bnd( 0, p )

	lin_solve(  0, p, div, 1, 4 )

	for i in range(1,N+1):
		for j in range(1,N+1):
			u[IX(i,j)] -= 0.5*N*(p[IX(i+1,j)]-p[IX(i-1,j)])
			v[IX(i,j)] -= 0.5*N*(p[IX(i,j+1)]-p[IX(i,j-1)])

	set_bnd( 1, u )
	set_bnd( 2, v )

def velocity_step( u, v, u0, v0 ):

	#u0[IX(33,33)]=5.-10.*random.random()
	v0[IX(int(N/2+1),int(N/2+1))]=5.

	add_source ( u, u0 )
	add_source ( v, v0 )
	#print_array(u)

	if diff>0:
		u0,u=SWAP(u0,u)
		diffuse(  1, u, u0 )
		v0,v=SWAP(v0,v)
		diffuse(  2, v, v0 )

	project(  u, v, u0, v0 )
	#print_array(u)

	u0,u=SWAP(u0,u)
	v0,v=SWAP(v0,v)

	advect(  1, u, u0, u0, v0 )
	advect(  2, v, v0, u0, v0 )
	project(  u, v, u0, v0 )

	return u,v,u0,v0
