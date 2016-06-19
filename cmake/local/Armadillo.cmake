# In order to use locally compiled Armadillo library, use the following 
# settings.
#
# If inlcude path to armadillo is $HOME/usr/local/include/armadillo, 
# un-comment the following line:
# list(APPEND CMAKE_PREFIX_PATH $ENV{HOME}/usr/local)
#
# If inlcude path to armadillo is $HOME/apps/armadillo/include/armadillo, 
# un-comment the following line:
list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/apps/armadillo")
