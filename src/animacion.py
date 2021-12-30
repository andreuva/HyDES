import numpy as np 
import cv2
import glob
from tqdm import tqdm


def animacion(name='animation', frames='*', path_out='', fps=15, size=(1200,1200)):
	
	fourcc = cv2.VideoWriter_fourcc(*'XVID')
	out = cv2.VideoWriter(f'{path_out}{name}.avi' ,fourcc ,fps ,size)

	lista = sorted(glob.glob(frames))

	for imag in tqdm(lista):
		img = cv2.imread(imag)
		print(f'writting image {imag} to video')
		img = cv2.resize(img,size) #redimensionar
		out.write(img)

	out.release()
