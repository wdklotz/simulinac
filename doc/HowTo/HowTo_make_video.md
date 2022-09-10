## The making of my tracking videos:
### Methode <u>imageJ</u>
* generate a sequence of frames in PNG format with SIMULINAC/tracker.py
* import all frames into imageJ and create AVI video
* convert the AVI video on-line to MP4 (https://www.online-convert.com/)
* watch MP4 video with windows media player (works also well on android)

### Methode <u>imagemagick</u>
* download imagemagick
* use imagemagick tool: convert -delay 25 -loop 1 *.png frames.gif
* use [online converter](https://anyconv.com/de/png-in-avi-konverter/) to convert frames.gif to avi

### Methode <u>ffmpeg</u>
* download ffmpeg
* ffmpeg -framerate 5 -i poincare_cut_%03d.png frames-09092022.avi
* ffmpeg has tons of options!
* fffmpeg has [python API](https://pypi.org/project/ffmpeg-python/)
* for %03d see [Image file demuxer](http://ffmpeg.org/ffmpeg-all.html#image2-1)