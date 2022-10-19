import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#import mysql.connector
import time 
import os 
from scipy.optimize import curve_fit
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import concurrent.futures
from pprint import pprint

def gauss(x, mu, sigma, A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)
def bimodal(x, mu1, sigma1, A1, mu2, sigma2, A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)
def save_multi_image(filename):
   pp = PdfPages(filename)
   fig_nums = plt.get_fignums()
   figs = [plt.figure(n) for n in fig_nums]
   for fig in figs:
      fig.savefig(pp, format='pdf')
   pp.close()
def CalcGain(img_CCD16):#Recibe 1 imagen DE 1 CCD
    data = (img_CCD16).flatten()
    data = np.delete(data, np.where(data == 0))
    a = list(range(-250,400))
    #fig = plt.figure()
    y,x,_ = plt.hist(data, bins = a)            
    try:
        x = (x[1:]+x[:-1])/2 
        expected = (0, .6, 400, 300, .8, 250)
        params, cov = curve_fit(bimodal, x, y, expected)
        gain = [str(params[3]-params[0])] 
        return gain
    except:
        print("Can't Fit Model")
        return [-1]
def GetSingleCCDImage(hdul,LTA_channel,ColInit,NCOL,tamxpimg,tamy,ccdncol):
    MuxedImage=hdul[LTA_channel].data
    step=112
    LastCol=ColInit+(NCOL-1)*step
    indexCol=list(range((ColInit),LastCol+step,step))
    print("Rango: ",indexCol[0],indexCol[1],indexCol[-1])
    DeMuxedImage=MuxedImage[:, indexCol]
    for p in range(tamy):
        Offset=np.mean(DeMuxedImage[p,(tamxpimg-int(NCOL-ccdncol/2)):tamxpimg])
        DeMuxedImage[p,:]=DeMuxedImage[p,:]-Offset
    return DeMuxedImage #return demuxed image

############################################################################################
if __name__ == '__main__': 
    directorio='/home/oscar/Documentos/Oscura/OscuraImagenes/ImagenesRepetidas/'
    directorio_guardado='/home/oscar/Documentos/Oscura/ImagenesProcesadas/'
    contenido = os.listdir(directorio)
    #Comienzo script
    st = time.time()
    img_CCD16=fits.HDUList([]) #Se crearan nro_imagenes*nsamp/16 (Ej: 5*112/6=35)
    gain_list=[]


    for item in range(len(contenido)):
        hdulist = fits.open(directorio+contenido[item])
        tamx=int(hdulist[4].header['NAXIS1']) #134400
        tamy=int(hdulist[4].header['NAXIS2']) #Varia 
        nsamp=int(hdulist[4].header['NSAMP']) #112
        ncol=int(hdulist[4].header['NCOL'])
        ccdncol=int(hdulist[4].header['CCDNCOL'])
        scidata = hdulist[4].data
        tamxpimg=int(tamx/nsamp)
        div=16 #Variable multiplo de 2^n 
        #img_parcial= np.zeros((tamy,tamxpimg),'i4')       #Se crean las imagenes parciales 
        gain=[]
        datos=[]
        contador=0
        for j in range(int(nsamp/16)):  #Se recorre nsamp/16=7 veces por cada imagen    
            for i in range(div):                               #Se recorre 16 veces para ir agarrando 16 CCD 
                img_parcial=GetSingleCCDImage(hdulist,4,i+div*j,ncol,tamxpimg,tamy,ccdncol)
                datos.append(img_parcial)
                img_CCD16.append(fits.ImageHDU(img_parcial))
            #Proceso de guardado
            directorio_demux=directorio_guardado+"Demuxed_"+contenido[item]+"/"
            if not os.path.exists(directorio_demux):
                os.makedirs(directorio_demux)
            nombre=str(directorio_demux+"MCM"+str(j+1)+"_Demuxed_"+contenido[item]+"_PROC.fits") 
            img_CCD16.writeto(nombre,overwrite=True)
            img_CCD16.clear()
        with concurrent.futures.ProcessPoolExecutor() as executor:
            result = executor.map(CalcGain, datos)
        print(tuple(result))           
        datos.clear()
        #gain.clear()    
        hdulist.close()

        #filename = directorio_demux=directorio_guardado+"Histogram112.pdf"
        #save_multi_image(filename)


    et = time.time()
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')