# simulinac
python simulator for linac structures.

I worked on a user friendly version of my code. Here it is!!!

Instructions:

1) you have to clone the code from github like that:

$git clone https://github.com/wdklotz/simulinac

then change directory:

$cd simulinac

then run the demo:

$python ./fodo.py

2) you must have python 3 (I tested with 3.4)

3) you need the modules: numpy, pylab, yaml (I use pyYaml)

4) there is a template input file fodo_template.yml. copy it and modify at your will

5) WARNING: you have to follow strictly the YAML syntax and don't remove (you can rename them) the label entries in the segments: and lattice: blocks. it will not parse correctly when you replace them. the label entries for segments are referenced in the lattice definition.

6) in the template I defined 2 identical segments for demonstration. you can define as many different segments as you like

7) in the lattice block you can have groups of segments in [...] brackets. the 1st integer in the list is mandatory and tells the parser to repeat the list behind as often as its value, i.e.  [10, SEG, SEGX] generates 10 times the sequence [SEG, SEGX]

8) the flags: parameters: and elements: blocks should be self explaining with the comments I added.

-- have fun --
