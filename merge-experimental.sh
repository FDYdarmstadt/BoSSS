#
# performes the splitting of the 'public' subtree from the fuul repository
# must be called from the public repo
#
git remote add root git@git.rwth-aachen.de:bosss1/experimental.git
git fetch root
git branch
git checkout -f -b staging-branch root/master
git subtree split -P public -b root-public
git checkout master
git merge -X theirs root-public
git branch -d root-public staging-branch
git remote remove root
git branch