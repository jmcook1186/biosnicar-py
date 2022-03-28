import React, { useState } from 'react';
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

    await fetch("http://localhost:5000/app/model", 
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

    <header style={{position:"absolute", top: 25, left: 100, color: 'black', fontSize: 80}}>
    <b>biosnicar</b>
    </header>

    <p style={{position:"absolute",left:100, top: 85, color: 'black', fontSize: 24}}>
    snow and ice albedo model</p>


    <p style={{position:"absolute",left:100, top: 180, color: 'black', fontSize: 22}}>
    Enter input values below</p>

    <li style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:100, top: 250, color: 'black', fontSize: 18}}>
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
    position:"absolute", left:100, top: 300, color: 'black', fontSize: 18}}>
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
    position:"absolute", left:100, top: 350, color: 'black', fontSize: 18}}>
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
    position:"absolute", left:100, top: 400, color: 'black', fontSize: 18}}>
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
    position:"absolute", left:100, top: 450, color: 'black', fontSize: 18}}>
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
    position:"absolute", left:100, top: 500, color: 'black', fontSize: 18}}>
      <p>Glacier algae conc (cells/mL)&nbsp;&nbsp;&nbsp;&nbsp;</p>
      
    <input 
      type="number"
      min="0"
      max="10000000"
      step="1"
      value={glacier_algae}
      placeholder="Enter algae concentration"
      onChange={e => setGA(e.target.value)} />
      </li>

    
    <li style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:100, top: 550, color: 'black', fontSize: 18}}>
      <p>Snow algae conc (cells/mL)&nbsp;&nbsp;&nbsp;&nbsp;</p>
      
    <input 
      type="number"
      min="0"
      max="10000000"
      step="1"
      value={snow_algae}
      placeholder="Enter algae concentration"
      onChange={e => setSA(e.target.value)} />
      </li>

    <li style={{display: 'flex', justifyContent:'center', alignItems:'center', color: 'black',
    position:"absolute", left:100, top: 600, color: 'black', fontSize: 18}}>
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
    position:"absolute", left:100, bottom: 170, fontSize: 28, width: 180, height: 100}} 
    onClick={runModel}><b>Submit</b></button>}&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;


    <form method="get" action={outFile}>
    <button type="submit" style={{position:"absolute", bottom:170, left:300, width: 180, 
    height: 100, fontSize: 18}}><b>Download albedo data</b></button>
    </form>

    <form method="submit" action="https://biosnicar-go-py.readthedocs.io/en/latest/">
    <button type="submit" style={{position:"absolute", bottom:170, left:500, width: 180, 
    height: 100, fontSize: 18}}><b>Documentation</b></button>
    </form>

    <p style={{position:'absolute', bottom:0, right:20, fontSize:18}}><i>Copyright: Joseph Cook (github.com/jmcook1186) 2022</i></p>

    <img src={figure} alt='albedo' style={{position:"absolute", left:950, top: 250, width: 800}} />


    </div>

  );
}

export default App;
