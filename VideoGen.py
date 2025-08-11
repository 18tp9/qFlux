#video generation from png function
import os
import cv2
from PIL import Image
import glob

def generate_video(image_folder,vidOutPath,fps):
    #os.chdir(vidOutPath)
    listdirect = sorted(os.listdir(image_folder))
    images = [img for img in listdirect if img.endswith(".jpg") or
             img.endswith(".jpeg") or img.endswith(".png")or img.endswith(".bmp")]
   
    fourcc = cv2.VideoWriter_fourcc(*'mp4v') 

    # Array images should only consider 
    # the image files ignoring others if any 
    
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    
    # setting the frame width, height width 
    # the width, height of first image 
    height, width, layers = frame.shape
    OutPath = vidOutPath
    video = cv2.VideoWriter(OutPath, fourcc,fps, (width, height))

    # Appending the images to the video one by one 
    for image in images:
        video.write(cv2.imread(os.path.join(image_folder,image)))
        
    # Deallocating memories taken for window creation 
    cv2.destroyAllWindows()
    video.release()  # releasing the video generated

def make_gif(gif_path,image_folder,dur):

    outPath = gif_path
    images = glob.glob(f"{image_folder}/*.png")
    images.sort()
    frames = [Image.open(image) for image in images]
    frame_one = frames[0]
    frame_one.save(outPath, format="GIF", append_images=frames[1:200:],save_all=True, duration=dur, loop=0)

def get_frames(filename):
    cap = cv2.VideoCapture(filename)
    ret, frame = cap.read()
    
    return frame

path_name = (r'C:/Users/pende/Documents/Summer22/Introduction/rhoAnimation')
# make_gif('../Figures/gifs/iter0.gif',path_name,1)
# make_gif('vids/dyePlots.gif','plots_dye',300)




# generate_video('stdA00Rough_sideByside','movies/stdA00R_sideByside.mp4',10)
# make_gif('movies/processed.gif','procTest',100)


# mamidir = ('Vvelo_depth')

# for kiddiedir in os.listdir(mamidir):

#     moviePath = os.path.join('movies',mamidir)#,kiddiedir)
    
#     try:
#         os.makedirs(moviePath)
#     except:
#         None
#     generate_video(os.path.join(mamidir,kiddiedir), 'movies/Vvelo_depth/'+ kiddiedir + '.mp4',10)
generate_video('FconvPNGs','Fdiff_Diff_kPhi_1_long.mp4',20)