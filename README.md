# simulinac
python simulator for linac structures.

I worked on a user friendly version of my code. Here it is!!!

Instructions:

1) You have to clone the code from github like that:

$git clone https://github.com/wdklotz/simulinac

then change directory:

$cd simulinac

then run the demo:

$python simu.py yml/ref_run.yml

2) You must have python 3 (I tested with 3.4).

3) You need the modules: numpy, pylab, yaml (I use pyYaml).

4) There is a reference input file yml/ref_run.yml. Copy it and modify at your will.

5) WARNING: You have to follow strictly the YAML syntax and don't remove (you can rename them) the label entries in the segments: and lattice: blocks. It will not parse correctly when you replace them. The label entries for segments are referenced in the lattice definition.

6) In the reference run I defined 2 Sections (sections: block) LE and HE. Only HE is used.  Element definitions (elements: block) have to be tagged by the section they belong to. You can define as many different segments as you like.

7) In the reference run I defined lattice segments (segments: block) SEG1H, SEG2H, SEG3H, RFGH, RFCAVH and MARK. You can define as many different segments from element combinations as you like.

8) In the lattice: block you can have groups of segments in [...] brackets. The 1st integer in the list is mandatory and tells the parser to repeat the list behind as often as its value, i.e. [10, SEG, SEGX] generates 10 times the sequence [SEG, SEGX].

9) The flags: parameters: and elements: blocks should be self explaining with the comments I added.

10) From the cloned repository read simu_manual.pdf which contains more information on how to prepare a lattice input.

11) You can get help from me by mail to wdklotz@gmail.com

-- have fun --
