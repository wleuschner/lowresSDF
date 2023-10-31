#!/usr/bin/env python

import time
import math
from sys import exit

import unicornhathd


print("""Unicorn HAT HD: Show a PNG image!

This basic example shows use of the Python Pillow library.

The tiny 16x16 bosses in lofi.png are from Oddball:
http://forums.tigsource.com/index.php?topic=8834.0

Licensed under Creative Commons Attribution-Noncommercial-Share Alike 3.0
Unported License.

Press Ctrl+C to exit!

""")

MAX_ITER = 255
MIN_DIST = 0.0
MAX_DIST = 100.0
EPSILON = 0.0001
WIDTH = 16
HEIGHT = 16

t = 0
torusToSphere = 0
startTime = time.time()

def lerp(a,b,alpha):
    alpha = max(min(alpha,1.0),0.0)
    return (1.0-alpha)*a+alpha*b

def dot(a,b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    
def absolute(a):
    return [math.fabs(a[0]),math.fabs(a[1]),math.fabs(a[2])]

def reflect(I,N):
    d = dot(N,I)
    dn2 = [2*d*N[0], 2*d*N[1], 2*d*N[2]]
    return [I[0] - dn2[0], I[1] - dn2[1], I[2] - dn2[2]]

def length(v):
    return math.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])

def normalize(v):
    l = length(v)
    return [v[0]/l,v[1]/l,v[2]/l]
    
def maxVec(a,b):
    return [max(a[0],b[0]),max(a[1],b[1]),max(a[2],b[2])]

def minVec(a,b):
    return [min(a[0],b[0]),min(a[1],b[1]),min(a[2],b[2])]

def sphereSDF(samplePoint, r):
    return length(samplePoint) - r
    
def torusSDF(samplePoint):
    q = [math.sqrt(samplePoint[0]*samplePoint[0]+samplePoint[1]*samplePoint[1])-0.06,samplePoint[2]]
    return math.sqrt(q[0]*q[0]+q[1]*q[1])-0.03
    
def cubeSDF(p, b):
    a = absolute(p)
    q = [a[0]-b[0],a[1]-b[1],a[2]-b[2]]
    c = maxVec(q,[0.0,0.0,0.0])
    d = min(max(q[0],max([q[1],q[2]])),0.0)
    return length(c) + d
    
def octahedronSDF(p, s):
    a = absolute(p)
    return (a[0]+a[1]+a[2]-s)*0.57735027
    
def sceneSDF(samplePoint):
    global startTime
    global torusToSphere
    
    angle = t
    samplePointRot = [samplePoint[0]*math.cos(angle) - samplePoint[2]*math.sin(angle), samplePoint[1], samplePoint[0] * math.sin(angle) + samplePoint[2] * math.cos(angle)]
    b = [0.1, 0.1,0.1]
    
    d1 = cubeSDF(samplePointRot,b)
    d2 = sphereSDF(samplePointRot,math.fabs(math.sin(t))*0.12)
    return max(d1,-d2)
    #p = samplePointRot
    #d1 = cubeSDF(p,b)
    #d2 = (math.sin(80*p[0]) * math.sin(80*p[1]) * math.sin(80*p[2])) * 0.01
    #return (d1+d2)
    #return octahedronSDF(samplePointRot, 0.075)
    #return cubeSDF(samplePointRot,[0.075, 0.075,0.1])
    
    morphTime = ((t-startTime)*1000*1000)/(2*1000*1000)
    ret = 0
    if torusToSphere == 0:
        ret = lerp(torusSDF(samplePointRot),sphereSDF(samplePointRot,0.1),morphTime)
    elif torusToSphere == 1:
        ret = lerp(sphereSDF(samplePointRot,0.1),cubeSDF(samplePointRot,b),morphTime)
    elif torusToSphere == 2:
        ret = lerp(cubeSDF(samplePointRot,b),torusSDF(samplePointRot),morphTime)
    
    if(morphTime>1.0):
        startTime = time.time()
        torusToSphere = (torusToSphere + 1) % 3
        print(torusToSphere)
    return ret
    
def shortesDistanceToSurface(eye,vdir,start,end):
    depth = start
    for i in range(MAX_ITER):
        offset = [eye[0]+depth*vdir[0],eye[1]+depth*vdir[1],eye[2]+depth*vdir[2]]
        dist = sceneSDF(offset)
        if dist<EPSILON:
            return depth
        depth = depth + dist
        if(depth >= end):
            return end
    return depth

def rayDirection(fov,size,fragCoord):
    xy = [fragCoord[0],fragCoord[1]]
    z = size[1] / math.tan(math.radians(fov) / 2.0)
    length = math.sqrt(xy[0]*xy[0]+xy[1]*xy[1]+z*z)
    return normalize([xy[0],xy[1],-z])
    
def estimateNormal(p):
    n = normalize([
    sceneSDF([p[0] + EPSILON, p[1], p[2]]) - sceneSDF([p[0] - EPSILON, p[1], p[2]]),
    sceneSDF([p[0], p[1] + EPSILON, p[2]]) - sceneSDF([p[0], p[1] - EPSILON, p[2]]),
    sceneSDF([p[0], p[1], p[2] + EPSILON]) - sceneSDF([p[0], p[1], p[2] - EPSILON])
    ])
    return n
    
def phongContribForLight(k_d, k_s, alpha, p, eye, lightPos, lightIntensity):
    N = estimateNormal(p)
    L = normalize([lightPos[0] - p[0],lightPos[1] - p[1],lightPos[2] - p[2]])
    V = normalize([eye[0] - p[0],eye[1] - p[1],eye[2] - p[2]])
    R = normalize(reflect([-L[0],-L[1],-L[2]],N))
    
    dotLN = dot(L, N)
    dotRV = dot(R, V)
    
    if dotLN < 0.0:
        return [0.0,0.0,0.0]
        
    if dotRV < 0.0:
        diff = [k_d[0] * dotLN,k_d[1] * dotLN,k_d[2] * dotLN]
        return [lightIntensity[0] * diff[0], lightIntensity[1] * diff[1], lightIntensity[2] * diff[2]]
    
    diffSpec = [k_d[0] * dotLN + k_s[0] * math.pow(dotRV, alpha),
                k_d[1] * dotLN + k_s[1] * math.pow(dotRV, alpha),
                k_d[2] * dotLN + k_s[2] * math.pow(dotRV, alpha)]
    return [lightIntensity[0] * diffSpec[0], lightIntensity[1] * diffSpec[1], lightIntensity[2] * diffSpec[2]]

def phongIllumination(k_a, k_d, k_s, alpha, p, eye):
    ambientLight = [0.5,0.5,0.5]
    color = [k_a[0] * ambientLight[0],k_a[1] * ambientLight[1],k_a[2] * ambientLight[2]]    
    
    light1Pos = [0.2 * math.sin(0.37 * t),0.2,0.2 * math.cos(0.37 * t)]
    light1Intensity = [0.4,0.4,0.4]
    light1Contrib = phongContribForLight(k_d, k_s, alpha, p, eye, light1Pos, light1Intensity);
    
    color = [color[0] + light1Contrib[0], color[1] + light1Contrib[1], color[2] + light1Contrib[2]]
    
    light2Pos = [0.2 * math.sin(0.37 * t), 0.2 * math.cos(0.37 * t), 0.2]
    light2Intensity = [0.4, 0.4, 0.4]
    light2Contrib = phongContribForLight(k_d, k_s, alpha, p, eye, light2Pos, light2Intensity)
    
    color = [color[0] + light2Contrib[0], color[1] + light2Contrib[1], color[2] + light2Contrib[2]]
    
    return color

unicornhathd.rotation(0)
unicornhathd.brightness(1.0)

width, height = unicornhathd.get_shape()


K_a = [0.2,0.2,0.2]
K_d = [0.7,0.2,0.2]
K_s = [1.0,1.0,1.0]
shininess = 10.0

try:
    while True:
        t = time.time()
        K_d[0] = math.fabs(math.sin(t))
        K_d[1] = math.fabs(math.cos(t))
        for y in range(HEIGHT):
            cy = (y - HEIGHT/2.0) / HEIGHT
            for x in range(WIDTH):
                cx = (x - WIDTH/2.0) / WIDTH
                vdir = rayDirection(45.0,[WIDTH,HEIGHT],[cx,cy])
                eye = [0.0,0.0,10.0]
                dist = shortesDistanceToSurface(eye,vdir,MIN_DIST,MAX_DIST)
                if dist > (MAX_DIST - EPSILON):
                    unicornhathd.set_pixel(x, y, 0, 0, 0)
                else:
                    p = [eye[0] + dist * vdir[0],eye[1] + dist * vdir[1],eye[2] + dist * vdir[2]]
                    n = phongIllumination(K_a,K_d,K_s,shininess,p,eye)
                    unicornhathd.set_pixel(x, y, n[0] * 255, n[1] * 255, n[2] * 255)
        unicornhathd.show()


except KeyboardInterrupt:
    unicornhathd.off()
