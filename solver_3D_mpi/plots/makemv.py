import ffmpeg as fp
stream = fp.input("./img%03d.png", framerate=10)
stream = fp.output(stream, "anime.mp4", pix_fmt='yuv420p')

fp.run(stream)