
#Run emacs on background
alias eb="emacs &"
alias seb="emacs --with-profile spacemacs &"
alias d='ssh discovery'
alias nu50='ssh nuweb50'
alias ph='ssh phong'
alias matlab='/Applications/MATLAB_R2019a.app/bin/matlab -nodesktop -nosplash'
alias rodeo='/Applications/Rodeo.app/Contents/MacOS/Rodeo'
alias PDFconcat='/System/Library/Automator/Combine\ PDF\ Pages.action/Contents/Resources/join.py'
alias jn='jupyter notebook'
alias jl='jupyter lab'
alias py36='source activate py36'
alias julia='exec /Applications/Julia-1.6.app/Contents/Resources/julia/bin/julia'
alias psh='pipenv shell'
alias j='julia'
alias jb='jupyter-book'
alias b='jb build ../Thesis'



#Add Path
export PATH="/Applications/Emacs.app/Contents/MacOS:$PATH"
export PATH="/Applications/Julia-1.6.app/Contents/Resources/julia/bin/:$PATH"


export GUROBI_HOME="/Library/gurobi902/mac64"

# Bash completion for jlpkg
# See https://github.com/fredrikekre/jlpkg for details.
if [[ -f ~/.bash_completion.d/jlpkg-completion.bash ]]; then
    . ~/.bash_completion.d/jlpkg-completion.bash
fi
