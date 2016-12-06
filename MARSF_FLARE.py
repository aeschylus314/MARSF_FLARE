import scipy as sp
from scipy.interpolate import InterpolatedUnivariateSpline
import scipy.interpolate
import matplotlib.pyplot as plt


def main():
  b0 = 0.474180489800000
  r0 = 0.858785390900000
  n = 3
  deltaphi_rad = 0.0 * (sp.pi/180)
  part = './' 
  bplas_u_path = part+'bplas_resist_resp_upper' # change here, 'vac' -> 'resist_resp'
  bplas_l_path = part+'bplas_resist_resp_lower' # change here, 'vac' -> 'resist_resp'
  rmzm_path = part+'rmzm_geom'
  profeq_path = part+'profeq'

  R_min = 0.1
  R_max = 1.9

  Z_min = -1.4
  Z_max = 1.4


  plot_or_save = 1 # change here. 0 = plot, 1 = save on rectangular grid
  if plot_or_save==0:
    nchi = 240
    rz = rzcoords(rmzm_path, nchi)
    jc = jacobian(rz)
    bplas_u = bplasma(bplas_u_path,rz,jc)
    bplas_l = bplasma(bplas_l_path,rz,jc)
    BN = (bplas_u.bn + bplas_l.bn*sp.exp(-1j*deltaphi_rad))*b0
    plt.contourf(rz.R*r0,rz.Z*r0, BN.imag, 120)
    plt.plot(rz.R[rz.Ns_plas,:]*r0,rz.Z[rz.Ns_plas,:]*r0,'-k', linewidth=2.0)
    plt.plot([R_min,R_min],[Z_max,Z_min],'-k',linewidth=2.0)
    plt.plot([R_max,R_max],[Z_max,Z_min],'-k',linewidth=2.0)
    plt.plot([R_min,R_max],[Z_min,Z_min],'-k',linewidth=2.0)
    plt.plot([R_min,R_max],[Z_max,Z_max],'-k',linewidth=2.0)
    plt.colorbar()
    plt.show()

  if plot_or_save==1:
    nchi= 240

    numR = 100
    numZ = 100
    numPhi = 60

    rz = rzcoords(rmzm_path, nchi)
    jc = jacobian(rz)
    bplas_u = bplasma(bplas_u_path,rz,jc)
    bplas_l = bplasma(bplas_l_path,rz,jc)
    
    BR = (bplas_u.Br + bplas_l.Br*sp.exp(-1j*deltaphi_rad))*b0
    BZ = (bplas_u.Bz + bplas_l.Bz*sp.exp(-1j*deltaphi_rad))*b0
    BP = (bplas_u.Bphi + bplas_l.Bphi*sp.exp(-1j*deltaphi_rad))*b0

    R1 = rz.R*r0
    Z1 = rz.Z*r0

    BR1 = BR
    BZ1 = BZ
    BP1 = BP

    R_rect = sp.linspace(R_min, R_max, numR)
    Z_rect = sp.linspace(Z_min, Z_max, numZ)

    R_grid, Z_grid = sp.meshgrid(R_rect,Z_rect)

    BR_rect = scipy.interpolate.griddata((R1.ravel(), Z1.ravel()), BR1.ravel(), (R_grid, Z_grid), method = 'linear')
    BR_rect = BR_rect.reshape((numZ,numR))

    BZ_rect = scipy.interpolate.griddata((R1.ravel(), Z1.ravel()), BZ1.ravel(), (R_grid, Z_grid), method = 'linear')
    BZ_rect = BZ_rect.reshape((numZ,numR))

    BP_rect = scipy.interpolate.griddata((R1.ravel(), Z1.ravel()), BP1.ravel(), (R_grid, Z_grid), method = 'linear')
    BP_rect = BP_rect.reshape((numZ,numR))
    #sp.savetxt('B_field_rectangular.txt', sp.transpose([R_grid.ravel(), Z_grid.ravel(), BR_rect.real.ravel(), BR_rect.imag.ravel(), BZ_rect.real.ravel(), BZ_rect.imag.ravel(), BP_rect.real.ravel(), BP_rect.imag.ravel()]), fmt='%.12e', header = "File format columns: R Z Re{Br} Im{Br} Re{Bz} Im{Bz} Re{Bp} Im{Bp} \nAn example plotting script is included in mogui_load_BRZ.py. RZ grid info is \n numR = %d, numZ = %d, Rmin = %0.4f, Rmax = %0.4f, Zmin = %0.4f, Zmax = %0.4f"%(numR,numZ, R_min, R_max, Z_min, Z_max))
    # Layout now defines the boxes for the outputs
    layout = sp.ndarray(shape=(numR, numZ, numPhi), dtype=complex)
    print(layout[0,0,0])
    print(BR_rect[0,0])
    BRoutcom = layout
    BZoutcom = layout
    BPoutcom = layout
    
    for x in xrange(0, numR):
        for y in xrange(0, numZ):
            for z in xrange(0, numPhi):
                BRoutcom[x,y,z]=BR_rect[y,x]*sp.exp(-1*1j*z*2*sp.pi/(numPhi*n))
    for x in xrange(0, numR):
        for y in xrange(0, numZ):
            for z in xrange(0, numPhi):
                BZoutcom[x,y,z]=BZ_rect[y,x]*sp.exp(-1*1j*z*2*sp.pi/(numPhi*n))
    for x in xrange(0, numR):
        for y in xrange(0, numZ):
            for z in xrange(0, numPhi):
                BPoutcom[x,y,z]=BP_rect[y,x]*sp.exp(-1*1j*z*2*sp.pi/(numPhi*n))
    
    BRoutreal=BRoutcom.real
    BZoutreal=BZoutcom.real
    BPoutreal=BPoutcom.real
    with file('MARSB_r.dat', 'w') as outfile1:
            outfile1.write('! MARS-F Output for BSpline input: Br \n')
            for data_slice in BRoutreal:
                sp.savetxt(outfile1, data_slice)
            outfile1.close()
    with file('MARSB_z.dat', 'w') as outfile2:
            outfile2.write('! MARS-F Output for BSpline input: Bz \n')
            for data_slice in BZoutreal:
                sp.savetxt(outfile2, data_slice)
            outfile2.close()
    with file('MARSB_phi.dat', 'w') as outfile3:
            outfile3.write('! MARS-F Output for BSpline input: Bphi \n')
            for data_slice in BPoutreal:
                sp.savetxt(outfile3, data_slice)
            outfile3.close()

class rzcoords():

  def __init__(self, path, nchi):
    rmzm = sp.loadtxt(path)
    Nm0 = rmzm[0,0] # Num. poloidal harmonics for equilibrium quantities (not necessarily same as for perturbation quantities, but should be).
    Ns_plas = rmzm[0,1] # Num. radial points in plasma
    Ns_vac = rmzm[0,2] # Num. radial points in vacuum
    Ns = Ns_plas + Ns_vac
    R0EXP = rmzm[0,3]
    B0EXP = rmzm[1,3]
    s = rmzm[1:Ns+1, 0]
    RM = rmzm[Ns+1:,0] + 1j*rmzm[Ns+1:,1]
    ZM = rmzm[Ns+1:,2] + 1j*rmzm[Ns+1:,3]
    RM = RM.reshape((Nm0, Ns))
    RM = sp.transpose(RM)
    ZM = ZM.reshape((Nm0, Ns))
    ZM = sp.transpose(ZM)
    RM[:,1:] = 2*RM[:,1:]
    ZM[:,1:] = 2*ZM[:,1:]

    m = sp.arange(0,Nm0,1)
    chi = sp.linspace(-sp.pi, sp.pi,nchi)
    expmchi = sp.exp(sp.tensordot(m,chi,0)*1j)
    R = sp.dot(RM[:,:Nm0],expmchi)
    Z = sp.dot(ZM[:,:Nm0],expmchi)
    
    self.R = sp.array(R.real) # R coordinates
    self.Z = sp.array(Z.real) # Z coordinates
    self.Nchi = nchi
    self.Nm0 = Nm0
    self.Ns_plas = Ns_plas  # number of s points in plasma
    self.Ns_vac = Ns_vac    # number of s points in vacuum
    self.Ns = Ns            # total number of s points
    self.R0EXP = R0EXP      # normalisation length
    self.B0EXP = B0EXP      # normalisation magnetic field
    self.m = m              # equilibrium poloidal harmonics
    self.chi=chi            # poloidal angle coordinate
    self.s=s                # radial coordinate = sqrt(psi_pol)


class jacobian():

  # Decide what stuff is needed later, and add a self. in front of it.
  def __init__(self, rz):
    if not isinstance(rz, rzcoords):
      print("Must pass in coordinate system of type plotting_base.rzcoords")
      return

    self.Ns = rz.Ns # Used in jacobian.plot(), so pass in from rzcoords

    self.dRds = sp.copy(rz.R) 
    self.dZds = sp.copy(rz.R)
    self.dRdchi = sp.copy(rz.R)
    self.dZdchi = sp.copy(rz.R)
    self.jacobian = sp.copy(rz.R)

  # this is for the vacuum region. these are overwritten for the plasma region
  # just having a number to denote the boundary index might be simpler in the future.
  # Vac_start variable should be all that's needed. II is way too complicated
    II_start = int(rz.Ns_plas)-1; II_end = len(rz.R[:,0])
    II2_start = int(rz.Ns_plas); II2_end = len(rz.R[:,0])

    s0 = sp.copy(rz.s[II_start:II_end]); R0 = sp.copy(rz.R[II_start:II_end, :])
    chi0=sp.squeeze(sp.copy(sp.array(rz.chi))); Z0 = sp.copy(rz.Z[II_start:II_end, :])

    hs = 0.5*(s0[1:] - s0[:-1]).min(); hs = min(hs,  2e-5)
    hchi = 0.5*(chi0[1:] - chi0[:-1]).min(); hchi = min(hchi,  1e-4)
    s1 = s0-hs; s2 = s0+hs
    chi1 = chi0 - hchi;  chi2 = chi0 + hchi

  # compute dR/ds using R(s,chi)
    R1 = sp.zeros(sp.shape(R0))
    R2 = sp.zeros(sp.shape(R0))
    for i in range(rz.Nchi):
      R1[:,i] = InterpolatedUnivariateSpline(s0,R0[:,i], bbox=[s1[0], s0[-1]])(s1)
      R2[:,i] = InterpolatedUnivariateSpline(s0,R0[:,i], bbox=[s0[0], s2[-1]])(s2)
    self.dRds[II_start:II_end,:] = (R2 - R1)/(2*hs)

  # compute dZ/ds using Z(s,chi) 
    Z1 = sp.zeros(sp.shape(Z0))
    Z2 = sp.zeros(sp.shape(Z0))
    for i in range(rz.Nchi):
      Z1[:,i] = InterpolatedUnivariateSpline(s0,Z0[:,i], bbox=[s1[0], s0[-1]])(s1)
      Z2[:,i] = InterpolatedUnivariateSpline(s0,Z0[:,i], bbox=[s0[0], s2[-1]])(s2)
    self.dZds[II_start:II_end,:] = (Z2 - Z1)/(2*hs)

  #  compute dR/dchi using R(s,chi) 
    R1 = sp.zeros(sp.shape(R0))
    R2 = sp.zeros(sp.shape(R0))
    for i in range(int(rz.Ns_vac)+1):
      R1[i,:] = InterpolatedUnivariateSpline(chi0,R0[i,:], bbox=[chi1[0], chi0[-1]])(chi1)
      R2[i,:] = InterpolatedUnivariateSpline(chi0,R0[i,:], bbox=[chi0[0], chi2[-1]])(chi2) 
      self.dRdchi[i+II_start,:] = (R2[i,:] - R1[i,:])/(2*hchi) 

  # compute dZ/dchi using Z(s,chi) 
    Z1 = sp.zeros(sp.shape(Z0))
    Z2 = sp.zeros(sp.shape(Z0))
    for i in range(int(rz.Ns_vac)+1):
      Z1[i,:] = InterpolatedUnivariateSpline(chi0,Z0[i,:], bbox=[chi1[0], chi0[-1]])(chi1)
      Z2[i,:] = InterpolatedUnivariateSpline(chi0,Z0[i,:], bbox=[chi0[0], chi2[-1]])(chi2)
    self.dZdchi[II_start:II_end,:] = (Z2 - Z1)/(2*hchi)

  # Now do same calculations for plasma region
    II_start=0; II_end = rz.Ns_plas;

    s0 = sp.copy(rz.s[II_start:II_end]); R0 = sp.copy(rz.R[II_start:II_end, :])
    chi0=sp.squeeze(sp.copy(sp.array(rz.chi))); Z0 = sp.copy(rz.Z[II_start:II_end, :])

    hs = 0.5*(s0[1:] - s0[:-1]).min(); hs = min(hs,  2e-5)
    hchi = 0.5*(chi0[1:] - chi0[:-1]).min(); hchi = min(hchi,  1e-4)
    s1 = s0-hs; s2 = s0+hs
    chi1 = chi0 - hchi;  chi2 = chi0 + hchi

  # compute dR/ds using R(s,chi) 
    R1 = sp.zeros(sp.shape(R0))
    R2 = sp.zeros(sp.shape(R0))
    for i in range(rz.Nchi):
      R1[:,i] = InterpolatedUnivariateSpline(s0,R0[:,i], bbox=[s1[0], s0[-1]])(s1)
      R2[:,i] = InterpolatedUnivariateSpline(s0,R0[:,i], bbox=[s0[0], s2[-1]])(s2)
    self.dRds[II_start:II_end,:] = (R2 - R1)/(2*hs)

  # compute dZ/ds using Z(s,chi) 
    Z1 = sp.zeros(sp.shape(Z0))
    Z2 = sp.zeros(sp.shape(Z0))
    for i in range(rz.Nchi):
      Z1[:,i] = InterpolatedUnivariateSpline(s0,Z0[:,i], bbox=[s1[0], s0[-1]])(s1)
      Z2[:,i] = InterpolatedUnivariateSpline(s0,Z0[:,i], bbox=[s0[0], s2[-1]])(s2)
      self.dZds[:II_end,i] = (Z2[:,i] - Z1[:,i])/(2*hs)

  #  compute dR/dchi using R(s,chi)
    R1 = sp.zeros(sp.shape(R0))
    R2 = sp.zeros(sp.shape(R0))
    for i in range(int(rz.Ns_plas)):
      R1[i,:] = InterpolatedUnivariateSpline(chi0,R0[i,:], bbox=[chi1[0], chi0[-1]])(chi1)
      R2[i,:] = InterpolatedUnivariateSpline(chi0,R0[i,:], bbox=[chi0[0], chi2[-1]])(chi2)
    self.dRdchi[II_start:II_end,:] = (R2 - R1)/(2*hchi)

  # compute dZ/dchi using Z(s,chi) 
    Z1 = sp.zeros(sp.shape(Z0))
    Z2 = sp.zeros(sp.shape(Z0))
    for i in range(int(rz.Ns_plas)):
      Z1[i,:] = InterpolatedUnivariateSpline(chi0,Z0[i,:], bbox=[chi1[0], chi0[-1]])(chi1)
      Z2[i,:] = InterpolatedUnivariateSpline(chi0,Z0[i,:], bbox=[chi0[0], chi2[-1]])(chi2)
    self.dZdchi[II_start:II_end,:] = (Z2 - Z1)/(2*hchi)

    G11 = sp.square(self.dRds) + sp.square(self.dZds)
    G12 = sp.multiply(self.dRds, self.dRdchi) + sp.multiply(self.dZds, self.dZdchi)
    G22 = sp.square(self.dRdchi) + sp.square(self.dZdchi)
    G22[0,:]=G22[1,:]
    G33 = sp.square(rz.R)

    # Metrics elements
    self.G11 = G11
    self.G12 = G12
    self.G22 = G22
    self.G33 = G33

    self.jacobian = (-self.dRdchi*self.dZds + self.dRds*self.dZdchi)*rz.R
    self.jacobian[0,:] = self.jacobian[1,:]



class bplasma():

  def __init__(self, path, rz, jc):
  
    self.path = path
    bplasma = sp.loadtxt(self.path)
    Nm1 = bplasma[0,0] # Number of perturbation poloidal harmonics (should be same as equilibrium harmonics)
    self.bm1 = bplasma[Nm1+1:, 0] + 1j*bplasma[Nm1+1:, 1]
    self.bm2 = bplasma[Nm1+1:, 2] + 1j*bplasma[Nm1+1:, 3]
    self.bm3 = bplasma[Nm1+1:, 4] + 1j*bplasma[Nm1+1:, 5]

    self.bm1 = self.bm1.reshape((Nm1, rz.Ns))
    self.bm2 = self.bm2.reshape((Nm1, rz.Ns))
    self.bm3 = self.bm3.reshape((Nm1, rz.Ns))

    for i in range(len(self.bm2[:,1])):
      self.bm2[i, 1:] = self.bm2[i, :-1]
      self.bm3[i, 1:] = self.bm3[i, :-1]

    m = sp.array(bplasma[1:Nm1+1,0])

    expmchi = sp.exp(sp.tensordot(m,rz.chi,0)*1j)

    self.b1 = sp.dot(self.bm1.T,expmchi)
    self.b2 = sp.dot(self.bm2.T,expmchi)
    self.b3 = sp.dot(self.bm3.T,expmchi)

    self.bn = self.b1/sp.sqrt(jc.G22*jc.G33)

    self.m=m 

    self.Br = sp.divide(sp.multiply(self.b1, jc.dRds) + sp.multiply(self.b2, jc.dRdchi), jc.jacobian)
    self.Bz = sp.divide(sp.multiply(self.b1, jc.dZds) + sp.multiply(self.b2,jc.dZdchi), jc.jacobian)
    self.Bphi = sp.divide(sp.multiply(self.b3, rz.R), jc.jacobian)

    self.Br[0,:] = self.Br[1,:]
    self.Bz[0,:] = self.Bz[1,:]
    self.Bphi[0:2,:] = self.Bphi[3,:]

    self.AbsB = sp.sqrt(sp.square(sp.absolute(self.Br)) + sp.square(sp.absolute(self.Bz)) + sp.square(sp.absolute(self.Bphi)))

    self.Nm1=Nm1


if __name__=='__main__':
  main()


