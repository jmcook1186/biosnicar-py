import React, { useState } from 'react';
import DownloadLink from 'react-download-link'
import './App.css';
import background from "./assets/background.jpg";
import figure from "./outputs/albedo.jpg"
import outFile from './outputs/albedo.csv';

function App() {
  
  const [layer_type, setLayerType] = useState()
  const [dz, setThickness] = useState()
  const [reff, setRadius] = useState()
  const [rho, setDensity] = useState()
  const [bc, setBC] = useState()
  const [glacier_algae, setGA] = useState()
  const [snow_algae, setSA] = useState()
  const [zenith, setZen] = useState()


  async function runModel(){
    console.log("inside func");
    console.log(layer_type, dz, reff, rho, bc, glacier_algae, snow_algae, zenith)

    await fetch("http://localhost:5000/model", 
    {'method': 'POST', 
    body: JSON.stringify({'lyr_typ': layer_type, 'dz': dz, 
      'r_eff': reff, 'rho': rho, 'bc': bc, 'glacier_algae': glacier_algae,
      'snow_algae': snow_algae, 'zenith': zenith}),
    headers: {'Content-Type':'application/json'}
    });

  }

  return (

    <div 
    className="App" 
    style={{ backgroundImage: 'url(' + background + ')',
    backgroundPosition: 'center',
    backgroundSize: 'cover',
    backgroundRepeat: 'no-repeat',
    width: '100vw',
    height: '100vh',
    }}>

    <header style={{position:"absolute", top: 25, left: 25, color: 'black', fontSize: 70}}>
    <b>biosnicar</b>
    </header>

    <p style={{position:"absolute",left:28, top: 75, color: 'black', fontSize: 22}}>
    snow and ice albedo model</p>

    <p style={{position:"absolute",right:55, top: 25, color: 'black', fontSize: 22}}>
    <a href="https://biosnicar-go-py.readthedocs.io/en/latest/"><b>Documentation</b></a>
    </p>


    <li style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:25, top: 150, color: 'black', fontSize: 18}}>
      <p>Layer Type&nbsp;&nbsp;&nbsp;&nbsp;</p>

    <input 
      type="number"
      min="0"
      max="1"
      value={layer_type}
      placeholder="Enter layer type (0 or 1)"
      onChange={e => setLayerType(e.target.value)} />
      </li>

    
    
    <li style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:25, top: 190, color: 'black', fontSize: 18}}>
      <p>Thickness (meters)&nbsp;&nbsp;&nbsp;&nbsp;</p>
      
    <input 
      type="number"
      min="0"
      max="50"
      step="0.01"
      value={dz}
      placeholder="Enter layer thickness"
      onChange={e => setThickness(e.target.value)} />
      </li>


    <li style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:25, top: 230, color: 'black', fontSize: 18}}>
      <p>Radius (microns)&nbsp;&nbsp;&nbsp;&nbsp;</p>
      
    <input 
      type="number"
      min="100"
      max="20000"
      step="100"
      value={reff}
      placeholder="Enter grain/bubble radius"
      onChange={e => setRadius(e.target.value)} />
      </li>


    <li style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:25, top: 270, color: 'black', fontSize: 18}}>
      <p>Density (kg/m3)&nbsp;&nbsp;&nbsp;&nbsp;</p>
      
    <input 
      type="number"
      min="200"
      max="900"
      step="50"
      value={rho}
      placeholder="Enter column density"
      onChange={e => setDensity(e.target.value)} />
      </li>



    <li style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:25, top: 310, color: 'black', fontSize: 18}}>
      <p>Black carbon conc (ppb)&nbsp;&nbsp;&nbsp;&nbsp;</p>
      
    <input 
      type="number"
      min="0"
      max="100000"
      step="1"
      value={bc}
      placeholder="Enter BC concentration"
      onChange={e => setBC(e.target.value)} />
      </li>



    <li style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:25, top: 350, color: 'black', fontSize: 18}}>
      <p>Glacier algae conc (cells/mL)&nbsp;&nbsp;&nbsp;&nbsp;</p>
      
    <input 
      type="number"
      min="0"
      max="100000"
      step="1"
      value={glacier_algae}
      placeholder="Enter algae concentration"
      onChange={e => setGA(e.target.value)} />
      </li>

    
    <li style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:25, top: 390, color: 'black', fontSize: 18}}>
      <p>Snow algae conc (cells/mL)&nbsp;&nbsp;&nbsp;&nbsp;</p>
      
    <input 
      type="number"
      min="0"
      max="100000"
      step="1"
      value={snow_algae}
      placeholder="Enter algae concentration"
      onChange={e => setSA(e.target.value)} />
      </li>

    <li style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:25, top: 430, color: 'black', fontSize: 18}}>
      <p>Solar zenith angle (degrees) &nbsp;&nbsp;&nbsp;&nbsp;</p>
      
    <input 
      type="number"
      min="0"
      max="100000"
      step="1"
      value={zenith}
      placeholder="Enter zenith angle"
      onChange={e => setZen(e.target.value)} />
      </li>


    {<button style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:50, bottom: 50, fontSize: 28, width: 150, height: 80}} 
    onClick={runModel}><b>Submit</b></button>}&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;


    <a href={outFile} download="albedo.csv"> Download Here </a>

    <img src={figure} alt='albedo' style={{position:"absolute", left:550, top: 200}} />


    </div>

  );
}

export default App;
