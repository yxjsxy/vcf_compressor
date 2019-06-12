import contextlib, sys
python3 = sys.version_info.major >= 3
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import os
import codecs
from scipy.linalg import svd
import subprocess
import arithmetic_coding


type = int(input('please input the type number(0 for compress, 1 for decompress):'))
filein = input('please input the filepath:')
num = input('enter the chr number:')

if type == 0:
	with codecs.open(filein, "r",encoding='utf-8', errors='ignore') as file:
		x = 0
		for line in file:
			datain = line.split()
			if datain[0] != num:
				x += 1
			else:
				fmt = datain[8].split(':')
				break
	path = os.getcwd()
	path_cpr = path + '/%s_compressed' %filein
	if not os.path.exists(path_cpr):
		os.mkdir(path_cpr)
	for i in fmt:
		if not os.path.exists(path_cpr+'/'+i):
			os.mkdir(path_cpr+'/'+i)
	print('starting analysing vcf file')
	print('We have the following annotations:')
	for i in fmt:
		print(i+'\n')
	indi = 0
	if len(fmt) > 1:
		with open(filein,"r") as file:
			gt = filein +'.GT'
			foutgt = open(gt,"w")
			y  = 0
			fout = [0]*(len(fmt)-1)
			for i in range(1,len(fmt)):
				fout[i-1] = open(path_cpr+'/'+fmt[i]+'/'+fmt[i],"w")
			print(len(fout))
			for line in file:
				y += 1
				if y<=x:
					foutgt.write(line)
				else:
					datain = line.split()
					if 'GT' in datain[8]:
						datain[8] = 'GT'
					for i in range(9,len(datain)):
						datain[i] = datain[i].split(':')
						a = datain[i][0]
							
						for j in range(len(fout)):
							if j + 1 < len(datain[i]):
								fout[j].write(datain[i][j+1]+'\t')
							else:
								#fout[j].write(datain[i][j+1]+'\t')
								fout[j].write(' '+'\t')
						if a:
							datain[i] = a
					dataout = "\t".join(datain)
					foutgt.write(dataout+'\n')
					for i in range(len(fout)):
						fout[i].write('\n')
			file.close()
			foutgt.close()
			for i in range(len(fout)):
				fout[i].close()
		print('Start compressing GT.')
		
		subprocess.Popen(["./gtc","compress", "-o",path_cpr+'/'+'GT'+'compressed',filein])
		
		print('GT compression finished.')
		if 'GL' in fmt:
			print('Start compressing GL.')
			
			subprocess.call(['gzip',path_cpr+'/GL/GL'])
			
			#data = pd.read_csv(path_cpr+'/GL/GL',sep='\t',index_col=False,header=None)
			
			print('GL compression finished.')
		if 'DP' in fmt:
			print('Start compressing DP.')
			data = pd.read_csv(path_cpr+'/DP/DP',sep='\t',index_col=False,header=None)
			print(data.iloc[:10,:10])
			data = data.replace('.',0)
			data = data.fillna(-1)
			data = np.array(data,dtype='float')
			print(np.amax(data))
			row,col = data.shape[0],data.shape[1]
			u, s, v = svd(data[:row//2,:], full_matrices=False)
			comp = 3
			pca = PCA(n_components=comp)
			pca.fit(data[:row//2,:])
			pc = pca.components_.astype(int)
        		#u = u[:,:comp]
			sigma = np.zeros((data[:row//2,:].shape[1],data[:row//2,:].shape[1]))
			sigma[:data[:row//2,:].shape[1],:data[:row//2,:].shape[1]] = np.diag(s)
			sigma = sigma[:,:comp]
			v = v[:comp,:]
			data1_r_o = np.dot(data[:row//2,:],pc.T).dot(pc)
			data_r_svd = np.dot(np.dot(u,sigma),v)
			#print(data1_r_o == data_r_svd,data1_r_o[1:5,1:5],data_r_svd[1:5,1:5])
			#print(pc[1:5,1:5], v[1:5,1:5])
			ar_1 = np.dot(data[:row//2,:],pc.T)
			ar_2 = np.dot(data[row//2:,:],pc.T)
			data1_r = ar_1.dot(pc)
			data2_r = ar_2.dot(pc)
        		#u_r, s, v = svd(data[row//2:,:],full_matrices=False)
			#print('svd time: '+str(time3 - time2))
			#u1_r = u1[:,:comp].astype(int)
			#u_r = u[:,:comp].astype(int)
			#s_r = s[:comp].astype(int)
			#v_r = v[:,:comp].astype(int)
			#print(u1.shape,u2.shape,u1_r.shape,u2_r.shape,s.shape,s_r.shape,v.shape,v_r.shape)
			a_r = np.vstack((ar_1,ar_1))
			#s_dia = np.diag(s_r)
			#data_r_svd = np.dot(u_r,s_dia).dot(v_r.T)
			data_r = np.vstack((data1_r,data2_r)).astype(int)
			error = data - data_r
			pdans = pd.DataFrame(error)
			print('error matrix formed.')
			outputfile = path_cpr+'/DP/compressed'
			os.mkdir(outputfile)
			np.savez_compressed(path_cpr+'/DP/matrix.npz',a=a_r,b=pc)
			#np.save(path_cpr+'/DP/PC',pc)
			#np.save(path_cpr+'/S',s_dia)
			print('U,V,S finished.')
			#pdans.to_csv(path_cpr+'/error',sep='\t',compression='gzip')
			#for i in ['A','PC','error']:
			#	with open(path_cpr+'/DP/'+i,'rb') as inp, contextlib.closing(arithmetic_coding.BitOutputStream(open(outputfile+'/'+i, "wb"))) as bitout:
			#		Encoder = arithmetic_coding.ArithmeticEncoder(32, bitout)
			#		Encoder.compress(inp, outputfile,bitout)
			np.savez_compressed(path_cpr+'/DP/error_matrix.npz',e=error)
			print('DP compression finished.')
	if len(fmt) == 1:
		print('Start compressing GT, which is the only annotation.')
		subprocess.Popen(["./gtc","compress", "-o",path_cpr+'/'+'GT'+'compressed',filein])
if type == 1:
	files = os.listdir(filein)
	findgt = 0
	for i in range(len(files)):
		path = os.path.join(filein,files[i])
		if os.path.isfile(path) and 'GT' in path and findgt==0:
			findgt = 1
			gt = path
			print('Start decompresseing GT.')
			subprocess.Popen(["./gtc","view","-o","chr"+num+".vcf",filein+"/GTcompressed"])
			print('GT decompression finished.')
	if findgt == 0:
		print('Error, no GT find.')
	findgl,finddp = 0,0
	for i in range(len(files)):
		path = os.path.join(filein,files[i])
		if os.path.isfile(path) and 'GL' in path:
			findgl = 1
			print('Start decompressing GL.')
			subprocess.call(['gunzip','-c',path+'/GL.gz','>','GL'])
			data_gl = np.load('GL')
			print('GL decompression fininshed.')
		if os.path.isfile(path) and 'DP' in path:
			finddp = 1
			matrix = np.load(path+'/matrix.npz')
			error = np.load(path+'/error_matrix.npz')
			a_r = matrix['a_r']
			pc = matrox['pc']
			data_r = np.dot(a_r,pc)
			data_dp = data_r + error
			print('DP decompression finished.')
	with open("chr"+num+".vcf","r") as file:
                x = 0
                for line in file:
                        datain = line.split()
                        if datain[0] != num:
                                x += 1
                        else:
                                fmt = datain[8].split(':')
	with open("chr"+num+".vcf","r") as file:
		fout = open('chr'+num+'.vcf.decompressed','w')
		y = 0 
		for line in file:
			y += 1
			if y<=x:
				fout.write(line+'\n')
			else:
				datain = line.split()
				if findgl == 1 :
					datain[8] += ':GL'
					for i in range(9,len(datain)):
						datain[i] += ':'+str(data_gl[y-x-1,i-9])
				if finddp == 1:
					datain[8] += ':DP'
					for i in range(9,len(datain)):
						datain[i] += ':'+str(data_dp[y-x-1,i-9])
				fout.write(datain.join('\t')+'\n')
	
