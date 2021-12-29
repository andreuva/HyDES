import numpy as np 
import cv2
import glob
from tqdm import tqdm

fourcc = cv2.VideoWriter_fourcc(*'XVID')
out = cv2.VideoWriter('animacion.avi' ,fourcc ,15 ,(1200,1200))

lista = sorted(glob.glob('*'))

for imag in tqdm(lista):
	img = cv2.imread(imag)
	# print(imag)
	# img = cv2.resize(img,(1000,1000)) #redimensionar
	out.write(img)

out.release()
