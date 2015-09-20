# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
fi

# User specific aliases and functions

export WM_SCHEDULER=ccache
export CCACHE_DIR=/vagrant/ccache4vm
