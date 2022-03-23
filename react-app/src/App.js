import React, { useState } from 'react';
import './App.css';
import background from "./assets/background.jpg";
import figure from "./outputs/albedo.jpg"

function App() {
  
  const [layer_type, setLayerType] = useState()
  const [dz, setThickness] = useState()
  const [reff, setRadius] = useState()
  const [rho, setDensity] = useState()

  // initial figure state
  //setFigure("./outputs/albedo_default.jpg")

  async function runModel(){
    console.log("inside func");
    console.log(layer_type, dz, reff, rho)

    await fetch("http://localhost:5000/model", 
    {'method': 'POST', 
    body: JSON.stringify({'lyr_typ': layer_type, 'dz': dz, 'r_eff': reff, 'rho': rho}),
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
    position:"absolute", left:25, top: 180, color: 'black', fontSize: 18}}>
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
    position:"absolute", left:25, top: 210, color: 'black', fontSize: 18}}>
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
    position:"absolute", left:25, top: 240, color: 'black', fontSize: 18}}>
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


    {<button style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:50, bottom: 50, fontSize: 28, width: 150, height: 80}} 
    onClick={runModel}><b>Submit</b></button>}&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;



    <img src={figure} alt='albedo' style={{position:"absolute", left:550, top: 200}} />


    </div>

  );
}

export default App;
