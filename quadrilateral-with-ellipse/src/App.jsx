import React from 'react';
import ReactDOM from 'react-dom/client';

const App = () => {
  return <h1>Hello, React!</h1>;
};

export default App;

const root = ReactDOM.createRoot(document.getElementById('root'));
root.render(
  <React.StrictMode>
    <App />
  </React.StrictMode>
);
