import cv2 as cv
import numpy as np
import glob
import tqdm
import os

# Create a VideoCapture object and read from input files
def create_video(path, filename, fps, size):
    fourcc = cv.VideoWriter_fourcc(*'XVID')
    out = cv.VideoWriter(os.path.join(path, filename), fourcc, fps, size)
    return out

# Write the frames to the video file
def write_video(out, frames):
    for i in tqdm.tqdm(range(len(frames))):
        out.write(frames[i])
    out.release()

# Read the frames from the image files
def read_frames(path):
    frames = []
    for filename in tqdm.tqdm(glob.glob(path)):
        frame = cv.imread(filename)
        frames.append(frame)
    return frames

# Create the video
def make_video(path, filename, fps, size):
    frames = read_frames(path)
    out = create_video(path, filename, fps, size)
    write_video(out, frames)

if __name__ == '__main__':
    path = 'results\\'
    filename = 'movie.avi'
    fps = 27
    size = (1920, 1080)
    make_video(path, filename, fps, size)