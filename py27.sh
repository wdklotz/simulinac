######### remove python3 from path
export PATH=$(echo "$PATH"|sed 's@/Users/klotz/anaconda3/bin:@@g')
echo'';echo 'PATH: '$PATH

######### PS1
export PS1="\[\033[01;32m\](py27)\[\033[00m\]"$PS1

