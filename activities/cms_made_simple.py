################################################################################
# CMS Made Simple
################################################################################

import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import mpl_toolkits.mplot3d.art3d as a3

import zipfile

import os

################################################################################
### I REALLY, REALLY HOPE THIS WORKS!!!!!!!!!!!!!!!
### 
### Yup! I think it works.
### Let's try this again....
###
### Third verse, same as the first....
###
################################################################################
def cms_made_simple_help(verbose=0):

    if verbose<1:
        print "Usage:\n"
        print "\tcollisions = get_collisions(f)"
        print "\nWhere f is a file (or zip file)"
        print "\n"
        print "\tjets,topjets,muons,electrons,met = collision"
        print "\tenergy,px,py,pz,csv = jet"
        print "\t#### EXPLANATION? AT LEAST OF CSV######"
        print "\tenergy,px,py,pz,nsub,minmass = topjet"
        print "\tenergy,px,py,pz = muon"
        print "\tenergy,px,py,pz = electron"
        print "\tpt,phi = met"
        print "\n"


    if verbose>=1:

        print "Usage:\n"
        print "\tcollisions = get_collisions(f)"
        print "\nWhere f is a file (or zip file)"
        print "\n"
        print "\tjets,topjets,muons,electrons,met = collision"
        print "\n"
        print "jets: Decay products of the quarks. Sprays of particles. (AK5Jets)"
        print "topjets: Jet from the top quark (CA8Jets)"
        print "muons & electrons- if a top decays leptonicly (the W decays to a lepton (muon or electron) and neutrino.\n"
        print "\tIt can also decay to a tau, but the tau cannot be detected because it decays too quickly)"
        print "met: Missing energy in the transverse plane. This is often used to identify if there was a neutrinos in the\n"
        print "\tevent, as they are not detected by CMS."
        print "\n"
        print "\tenergy,px,py,pz,csv = jet"
        print "\tenergy,px,py,pz,nsub,minmass = topjet"
        print "\tenergy,px,py,pz = muon"
        print "\tenergy,px,py,pz = electron"
        print "\tpt,phi = met"
        print "\n"
        print "jet:\n"
        print "\tenergy-energy of the particle/jet"
        print "\tpx-momentum in the x direction"
        print "\tpy-momentum in the y direction"
        print "\tpz-momentum in the z direction"
        print "\tcsv: Combined Secondary Vertex. This variable ranges from 0-1 (CHECK THIS!) and attempts to distinguish"
        print "\t\tbetween b-quarks which live longer and fly further before decaying, and lighter quarks which decay very"
        print "\t\tquickly. The closer this value is to 1, the greater confidence we have that the jet came from a b-quark."
        print "topjet:\n"
        print "\tenergy-energy of the particle/jet"
        print "\tpx-momentum in the x direction"
        print "\tpy-momentum in the y direction"
        print "\tpz-momentum in the z direction"
        print "\tnsub: Number of sub-jet structure. Depending on how boosted the top is, the algorithm can still look for"
        print "\t\tsubstructure which would correspond to the jets from the decay of the b-quark or W boson."
        print "\tminmass-"
        print "pt-momentum in the transverse plane.  Momentum that is orthongonal to the beam."
        print "phi-called the azimuthal angle.  It is the angle between the pt and x-axis if the beams are on the z. "
        print "\n"

################################################################################
################################################################################
def pretty_print(collision):

    jets,topjets,muons,electrons,photons,met = collision

    print "------- jets"
    for p in jets:
        energy,px,py,pz,csv = p
        print "energy:%8.5f px:%12.5f py:%12.5f pz:%12.5f csv:%12.5f" % (energy,px,py,pz,csv)
    print "------- top jets"
    for p in topjets:
        energy,px,py,pz,nsub,minmass = p
        print "energy:%8.5f px:%12.5f py:%12.5f pz:%12.5f nsub:%12.5f min mass:%12.5f" % (energy,px,py,pz,nsub,minmass)
    print "------- muons"
    for p in muons:
        energy,px,py,pz,q = p
        print "energy:%8.5f px:%12.5f py:%12.5f pz:%12.5f charge: %d" % (energy,px,py,pz,q)
    print "------- electrons"
    for p in electrons:
        energy,px,py,pz,q = p
        print "energy:%8.5f px:%12.5f py:%12.5f pz:%12.5f charge: %d" % (energy,px,py,pz,q)
    print "------- photons"
    for p in photons:
        energy,px,py,pz = p
        print "energy:%8.5f px:%12.5f py:%12.5f pz:%12.5f" % (energy,px,py,pz)
    print "------- met"
    for p in met:
        pt,phi = p
        print "pt:%8.5f phi:%8.5f" % (pt,phi)


################################################################################
def get_collisions(infilename,verbose=False):

    infile = None
    if zipfile.is_zipfile(infilename) is True:
        z = zipfile.ZipFile(infilename,'r')
        infile = z.open(z.namelist()[0],'r')
    else:
        infile = open(infilename)

    collisions = get_collisions_from_file_object(infile,verbose)

    return collisions


################################################################################
def ptetaphi_to_xyz(pt,eta,phi):
    px = pt*np.cos(phi)
    py = pt*np.sin(phi)
    pz = pt*np.sinh(eta)

    return px,py,pz


################################################################################
def get_collisions_from_file_object(infile,verbose=False):

    collisions = []

    not_at_end = True
    collision_count = 0
    new_collision = True
    while ( not_at_end ):

        ############################################################################
        # Read in one collision
        ############################################################################
        line = infile.readline()

        if collision_count%1000==0 and verbose:
            print "collision count: ",collision_count

        if line=="":
            not_at_end = False

        if line.find("Event")>=0:
            new_collision = True

        if new_collision==True:

            # Read in the jet info for this collision.
            jets = []
            line = infile.readline()
            njets = int(line)
            for i in xrange(njets):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                bquark_jet_tag = float(vals[4])
                jets.append([e,px,py,pz,bquark_jet_tag])

            # Read in the top jet info for this collision.
            topjets = []
            line = infile.readline()
            ntopjets = int(line)
            for i in xrange(ntopjets):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                nsub = float(vals[4])
                minmass = float(vals[5])
                topjets.append([e,px,py,pz,nsub,minmass])

            # Read in the muon info for this collision.
            muons = []
            line = infile.readline()
            nmuons = int(line)
            num_mu=0
            for i in xrange(nmuons):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                charge = int(float(vals[4]))
                muons.append([e,px,py,pz,charge])
                #muons.append([e,px,py,pz])
                num_mu+=1
                

            # Read in the electron info for this collision.
            electrons = []
            line = infile.readline()
            nelectrons = int(line)
            for i in xrange(nelectrons):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                charge = int(float(vals[4]))
                electrons.append([e,px,py,pz,charge])
                #electrons.append([e,px,py,pz])

            # Read in the photon info for this collision.
            #'''
            photons = []
            line = infile.readline()
            nphotons = int(line)
            for i in xrange(nphotons):
                line = infile.readline()
                vals = line.split()
                e = float(vals[0])
                px = float(vals[1])
                py = float(vals[2])
                pz = float(vals[3])
                photons.append([e,px,py,pz])
            #'''


            # Read in the information about the missing transverse energy (MET) in the collision.
            # This is really the x and y direction for the missing momentum.
            met = []
            line = infile.readline()
            nmet = int(line)
            for i in xrange(nmet):
                line = infile.readline()
                vals = line.split()
                #met_px = float(vals[0])
                #met_py = float(vals[1])
                met_pt = float(vals[0])
                met_phi = float(vals[1])
                met.append([met_pt,met_phi])

            new_collision = False
            collision_count += 1

            collisions.append([jets,topjets,muons,electrons,photons,met])

    return collisions

###############################################################################
'''
def get_collisions_from_zipped_file(infile,verbose=False):
    return get_collisions(infile)
'''


################################################################################
################################################################################
def draw_jet(origin=(0,0),angle=90,length=0.5,opening_angle=20,ntracks=5,show_tracks=False):

    lines = []
    patches = []

    # Edges of cone
    width_at_top = length*np.deg2rad(opening_angle)
    for side in [-1,1]:
        theta0 = np.deg2rad(angle+(side*opening_angle/2.0)) 
        x1 = length*np.cos(theta0)
        y1 = length*np.sin(theta0)
        #print x1,y1
        line = mlines.Line2D((origin[0],x1), (origin[1],y1), lw=2., alpha=0.4,color='red',markeredgecolor='red')
        lines.append(line)

    # End of cone
    arad = np.deg2rad(angle)
    center = (origin[0]+np.cos(arad)*length,origin[1]+np.sin(arad)*length)
    #print center
    p = mpatches.Ellipse(center, width_at_top+0.01, width_at_top/2.0,facecolor='red',alpha=0.4,edgecolor='gray',angle=abs(angle+90))
    patches.append(p)

    return patches,lines


    

################################################################################
################################################################################
def draw_jets(origins=[(0,0)],angles=[90],lengths=[0.5],opening_angles=[20],ntrackss=[5],show_trackss=[False]):

    alllines = []
    allpatches = []

    # Edges of cone
    for origin,angle,length,opening_angle,ntracks,show_tracks in zip(origins,angles,lengths,opening_angles,ntrackss,show_trackss):
        patches,lines = draw_jet(origin=origin,angle=angle,length=length,opening_angle=opening_angle,ntracks=ntracks,show_tracks=show_tracks)
        allpatches += patches
        alllines += lines


    return allpatches,alllines


    
################################################################################
################################################################################
def draw_line3D(origin=[(0,0,0)],pmom=[(1,1,1)],color='red',ls='-',lw=2.0,alpha=1.0):

    lines = []

    #print pmom
    for o,p in zip(origin,pmom):
        #x1 = p[0]
        #y1 = p[1]
        #z1 = p[2]
        x1 = p[2]
        y1 = p[0]
        z1 = p[1]
        #print x1,y1,z1
        line = a3.Line3D((o[0],x1),(o[1],y1),(o[0],z1), lw=lw, ls=ls, alpha=alpha,color=color,markeredgecolor=color)
        lines.append(line)

    return lines


################################################################################
################################################################################
def draw_beams():

    lines = draw_line3D(origin=[(0,0,-0.1),(0,0,0.1)],pmom=[(0,0,-200.0),(0,0,200.0)],color='red',lw=1,ls='--')

    return lines

################################################################################
################################################################################
def draw_jet3D(origin=[(0,0,0)],pmom=[(1,1,1)]):

    neworg = origin.copy()
    newmom = pmom.copy()

    offset = [[0.05,0.05,0.05],
              [0.05,0.05,-0.05],
              [0.05,-0.05,0.05],
              [0.05,-0.05,-0.05],
              [-0.05,0.05,0.05],
              [-0.05,0.05,-0.05],
              [-0.05,-0.05,0.05],
              [-0.05,-0.05,-0.05],
            ]

    offset = np.array(offset)
    offset *= 50

    for p in pmom:
        for o in offset:
            #print p.copy(),o
            pnew = p.copy() + o
            #print pnew
            newmom = np.vstack((newmom,pnew))
            neworg = np.vstack((neworg,(0,0,0)))

    lines = draw_line3D(origin=neworg,pmom=newmom,color='orange',lw=1)

    return lines

################################################################################
################################################################################
def draw_muon3D(origin=[(0,0,0)],pmom=[(1,1,1)]):

    lines = draw_line3D(origin=origin,pmom=pmom,color='blue',lw=5)

    return lines

################################################################################
################################################################################
def draw_electron3D(origin=[(0,0,0)],pmom=[(1,1,1)]):

    lines = draw_line3D(origin=origin,pmom=pmom,color='green',lw=2)

    return lines


################################################################################
################################################################################
def draw_photon3D(origin=[(0,0,0)],pmom=[(1,1,1)]):

    lines = draw_line3D(origin=origin,pmom=pmom,color='yellow',ls='-',lw=5)

    return lines

################################################################################
def draw_MET3D(origin=[(0,0,0)],pmom=[(1,1,1)]):

    lines = draw_line3D(origin=origin,pmom=pmom,color='red',ls='-',lw=9,alpha=0.5)

    return lines

################################################################################
################################################################################
################################################################################
def display_collision3D(collision,fig=None,ax=None):

    jets,topjets,muons,electrons,photons,met = collision

    lines = draw_beams()

    pmom = np.array(jets).transpose()[1:4].transpose()
    origin = np.zeros((len(jets),3))
    lines += draw_jet3D(origin=origin,pmom=pmom)

    pmom = np.array(muons).transpose()[1:4].transpose()
    origin = np.zeros((len(muons),3))
    lines += draw_muon3D(origin=origin,pmom=pmom)

    pmom = np.array(photons).transpose()[1:4].transpose()
    origin = np.zeros((len(photons),3))
    lines += draw_photon3D(origin=origin,pmom=pmom)

    pmom = np.array(electrons).transpose()[1:4].transpose()
    origin = np.zeros((len(electrons),3))
    lines += draw_electron3D(origin=origin,pmom=pmom)

    pmom = [ptetaphi_to_xyz(met[0][0],0.0,met[0][1])]
    origin = np.zeros((1,3))
    lines += draw_MET3D(origin=origin,pmom=pmom)

    if ax is None and fig is not None:
        ax = fig.add_subplot(1,1,1)
        ax = fig.gca(projection='3d')
        plt.subplots_adjust(top=0.98,bottom=0.02,right=0.98,left=0.02)
    
    if ax is not None:
        ax.cla() # Clear the axes

    if fig is None:
        fig = plt.figure(figsize=(7,5),dpi=100)
        if ax is None:
            ax = fig.add_subplot(1,1,1)
            ax = fig.gca(projection='3d')
            plt.subplots_adjust(top=0.98,bottom=0.02,right=0.98,left=0.02)

    for l in lines:
        ax.add_line(l)

    ax.set_xlim(-200,200)
    ax.set_ylim(-200,200)
    ax.set_zlim(-200,200)

    return lines,fig,ax

################################################################################
# Write function
################################################################################
def write_to_file(collisions,filename='default_collisions_output',do_zip=True):

    outfilename = "%s.dat" % (filename)
    outfile = open(outfilename,'w')

    #print "\tenergy,px,py,pz,csv = jet"
    #print "\tenergy,px,py,pz,nsub,minmass = topjet"
    #print "\tenergy,px,py,pz = muon"
    #print "\tenergy,px,py,pz = electron"
    #print "\tpt,phi = met"

    i = 0
    for collision in collisions:

        #print collision
        jets,topjets,muons,electrons,photons,met = collision
        
        if i%1000==0:
            print "Event: %d" % (i)

        output = "Event: %d\n" % (i)

        ############################################################################
        # Print out the not-top jets
        ############################################################################
        output += "%d\n" % (len(jets))
        for jet in jets:
            energy,px,py,pz,btag = jet
            output += "%-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n" % (energy,px,py,pz,btag)

        ############################################################################
        # Print out the top jets
        ############################################################################
        output += "%d\n" % (len(topjets))
        for jet in topjets:
            energy,px,py,pz,nsub,minmass = jet
            output += "%-10.4f %-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n" % (energy,px,py,pz,nsub,minmass)

        ############################################################################
        # Print out the muons
        ############################################################################
        output += "%d\n" % (len(muons))
        for muon in muons:
            energy,px,py,pz,q = muon
            output += "%-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n" % (energy,px,py,pz,q)

        ############################################################################
        # Print out the electrons
        ############################################################################
        output += "%d\n" % (len(electrons))
        for electron in electrons:
            energy,px,py,pz,q = electron
            output += "%-10.4f %-10.4f %-10.4f %-10.4f %-10.4f\n" % (energy,px,py,pz,q)

        ############################################################################
        # Print out the photons
        ############################################################################
        output += "%d\n" % (len(photons))
        for photon in photons:
            energy,px,py,pz = photon
            output += "%-10.4f %-10.4f %-10.4f %-10.4f\n" % (energy,px,py,pz)


        ############################################################################
        # Print out the met
        ############################################################################
        output += "%d\n" % (len(met))
        for m in met:
            pt,phi = m
            output += "%-10.4f %-10.4f\n" % (pt,phi)

        outfile.write(output)

        i += 1

    outfile.close()

    if do_zip:
        zipfilename = "%s.zip" % (filename)
        zf = zipfile.ZipFile(zipfilename,'w')
        zf.write(outfilename,compress_type=zipfile.ZIP_DEFLATED)
        zf.close()

        os.remove(outfilename)

    return 0
