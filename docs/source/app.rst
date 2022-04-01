*************
BioSNICAR App
*************

The BioSNICAR app is designed to be as easy as possible to use for non-programmers. Apart from running the app, absolutely no coding is required to generate albedo data
directly in a normal web browser. This is achieved using a Flask back-end that interacts with the BioSNICAR model code. This is a Flask server running on port 5000 that accepts user data in the form
of http POST requests sent from the front end. This data is used to run BioSNICAR and update an albedo figure and a csv file containing albedo data which are served by the front end.
The front end itself is a single page React app running on port 3000. The app has not yet been deployed to a public server, but it can be run locally after downloading the BioSNCAR source code
from the Github repository.

It is assumed that the user has a Python environment (venv or conda) that contains the basic BiOSNICAR dependencies as defined in `requirements.txt`.

Using the App
-------------

Intrustions for using the app are provided below, and there is also a walkthrough video provided in the project README. Please note that so far the app has only been tested on Linux (Ubuntu 20.04) and Firefox.

Installing node dependencies
----------------------------

The BioSNICAR app is composed of a Flask backend connected to a React frontend. The backend runs on the dependencies already installed in the Python environment described in the previous section, but the frontend requires some javascript packaged to be installed too. First check if nNdejs is already installed on your machien by opening a terminal and running:

`node -v`

If you already have Nodejs installed you will see a version number in the terminal. If not, download Nodejs for your operating system [here](https://nodejs.org/en/). We then also want to install yarn. Your Nodejs installation includes Node Package Manager (npm) which can now be used to install yarn. Open a terminal and run:

`npm install --global yarn`

Yarn can now be used to install the specific javascript packages needed to run the frontend. To do this navigate to the `app/src` directory and run:

`yarn`

You will see notifications for quite a few package installations. Now you have everything required to run the BioSNICAR app.

Run the app
-----------

The code for the Flask backend is in `~/app/api/app.py`. The code for the React front-end is in `app/src/`. 

In a terminal, navigate to the top-level BioSNICAR directory and run:


`python ./app/api/app.py`

This starts the Flask server running on `localhost:5000`. You can now ignore this and open a new terminal window (leaving the original terminal running in the background).

In the new terminal navigate to `/app/src/` and run:

`yarn start`

The BioSNICAR app will automatically open in your web browser. If it doesn't, or you accidentally close it,. you can access it again by navigating the browser to `localhost:3000`.


Get albedo data
----------------


Simply update the values in the input fields and press `Submit`. The spectral albedo plot and the broadband albedo value will update on the screen. You can download this data to a csv file by clicking `download data`. 

Video walkthrough
-----------------

The walkthrough video for the app is available at https://youtu.be/9wHMzZAB_do.