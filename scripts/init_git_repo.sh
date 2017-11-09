# Create Bitbucket repo, add scripts folder to versioning

curl --user mgloud@stanford.edu:<pass> https://api.bitbucket.org/1.0/repositories/ --data name=motrpac_rats --data is_private='true'
git init /srv/gsfs0/projects/montgomery/mgloud/projects/motrpac/rats 
git remote add origin https://mgloud@bitbucket.org/mgloud/motrpac_rats.git
git add *

