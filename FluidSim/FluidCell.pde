
final int N =128;
final int SCALE =10;
float t = 0;
int iter = 4;


class FluidCell {
  int size;
  float diffusion;
  float viscosity;
  float time;
  
  float[] Vx;
  float[] Vy;
  float[] density;
  
  float[] prevVx;
  float[] prevVy;
  float[] prevDensity;
  
  
  
  public FluidCell(float time, float diffusion, float viscosity) {
    this.size = N;
    this.diffusion = diffusion;
    this.viscosity = viscosity;
    this.time = time;
    
    this.Vx = new float[N*N];
    this.Vy =  new float[N*N];
    this.density =  new float[N*N];
    
    this.prevVx =  new float[N*N];
    this.prevVy = new float[N*N] ;
    this.prevDensity =  new float[N*N];
  }
  
  
  void addDye(int x, int y, float amt) {
    //add dye (impacts density) to specific spot in fluid cell
    this.density[index(x,y)] += amt;  
  }
  
  
  void addVelocity(int x, int y, float amtX, float amtY) {
    this.Vx[index(x,y)] += amtX;
    this.Vy[index(x,y)] += amtY;    
  }
  
  
  void diffuse(int b, float[] x, float[] prevX, float diffusion, float time) {
   //essentially solves navier-stokes eqn; diffuses dye into fluid
    float a = time * diffusion * (N-2)*(N-2);
    linSolve(b,x,prevX,a,1 + 4 * a);
  }
   
  
  void linSolve(int b, float[] x, float[] prevX, float a, float c) {
    //runs through array and sets each individual point in fluid to a combination of its neighbors, also resets bounds at each iter
    float cRecip = 1/c;
    for (int k = 0; k < iter; k++) { //track iters
      for (int i = 1; i < N-1; i++) { //track x coord
        for (int j = 1; j < N-1; j++) { //track y coord
          x[index(i, j)] =
          (prevX[index(i, j)]
          + a*(    x[index(i+1, j)]
          +x[index(i-1, j)]
          +x[index(i, j+1)]
          +x[index(i, j-1)]
          )) * cRecip;
        }
        
      }
      setBounds(b,x);
    }    
  }
  
  
  void project(float[] Vx, float[] Vy, float[] p, float[] div) {
    //net inflow of fluid = net outflow; maintain equilibrium
    for (int j = 1; j < N-1; j++) {
      for (int i = 1; i < N-1; i++) {
        div[index(i, j)] = -0.5f*(
                         Vx[index(i+1, j)]
                        -Vx[index(i-1, j)]
                        +Vy[index(i , j+1 )]
                        -Vy[index(i , j-1)]
                    )/N;
                p[index(i, j)] = 0;
      }
    }    
    setBounds(0,div);
    setBounds(0,p);
    linSolve(0,p,div,1,6);
    
    for (int j = 1; j < N-1; j++) {
      for (int i = 1; i < N-1; i++) {
      Vx[index(i, j)] -= 0.5f * (  p[index(i+1, j)]
                                                -p[index(i-1, j)]) * N;
                Vy[index(i, j)] -= 0.5f * (  p[index(i, j+1)]
                                                -p[index(i, j-1)]) * N;
               
      }    
    }    
    setBounds(1, Vx);
    setBounds(2, Vy);  
  }
  
  
  void advect(int b, float[] d, float[] prevD, float[] Vx, float[] Vy, float dt) { //WRONG WRONG WRONG WRONG 
  //moves stuff around! looks at each point, gets velocity, follows velocity backwards and sees where it lands. 
  // then takes weighted avg of cells around landing spot, and applies value to current cell
  //note that advection applies to the dye and the velocities themselves
    float i0, i1, j0, j1;
    
    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);
    
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    
    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;
  
    for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
       for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
          tmp1 = dtx * Vx[index(i, j)];
          tmp2 = dty * Vy[index(i, j)];
                
          x    = ifloat - tmp1; 
          y    = jfloat - tmp2;
                
                
          if(x < 0.5f) x = 0.5f; 
          if(x > Nfloat + 0.5f) x = Nfloat + 0.5f; 
          i0 = floor(x); 
          i1 = i0 + 1.0f;
          if(y < 0.5f) y = 0.5f; 
          if(y > Nfloat + 0.5f) y = Nfloat + 0.5f; 
          j0 = floor(y);
          j1 = j0 + 1.0f;
                
          s1 = x - i0; 
          s0 = 1.0f - s1; 
          t1 = y - j0; 
          t0 = 1.0f - t1;
                
          int i0i = (int) i0;
          int i1i = (int)i1;
          int j0i =(int) j0;
          int j1i = (int)j1;

                
          d[index(i, j)] = s0 * (t0 * prevD[index(i0i, j0i)] + t1 * prevD[index(i0i, j1i)]) +
                           s1 * (t0 * prevD[index(i1i, j0i)] + t1 * prevD[index(i1i, j1i)]);
          }
            
      }
      setBounds(b,d);
  
    }
  
  void setBounds(int b, float[] x) {
    //walls
   for(int i = 1; i < N - 1; i++) {
        x[index(i, 0  )] = b == 3 ? -x[index(i, 1  )] : x[index(i, 1  )];
        x[index(i,N-1)] = b == 3 ? -x[index(i, N-2)] : x[index(i, N-2)];
    }
  
    for(int j = 1; j < N - 1; j++) {
        x[index(0  , j)] = b == 1 ? -x[index(1  , j)] : x[index(1  , j)];
        x[index(N-1, j)] = b == 1 ? -x[index(N-2, j)] : x[index(N-2, j)];
    }

    
    //corners
    x[index(0,0)] = 0.5f * (x[index(1,0)] + x[index(0,1)]);
    x[index(0,N-1)] = 0.5f * (x[index(1,N-1)] + x[index(0,N-2)]);
    x[index(N-1,0)] = 0.5f * (x[index(N-2,0)] + x[index(N-1,1)]);
    x[index(N-1,N-1)] = 0.5f * (x[index(N-2,N-1)] + x[index(N-1,N-2)]); 
  }
  
  

  int index(int x, int y) {
    //get 2D location from 1D array
    x = constrain(x,0,N-1);
    y = constrain(y,0,N-1);
    return x + y * N;
  }
  
  
  void timeStep() {
    //diffuse x and y velocities based on viscosity and time
   diffuse(1, this.prevVx, this.Vx, this.viscosity, this.time);
   diffuse(2, this.prevVy, this.Vy, this.viscosity, this.time);

   
   //clean up to preserve equilibrium of non compressible fluid
   project(this.prevVx, this.prevVy, this.Vx, this.Vy);

  // advection on x and y velocities
   advect(1, this.Vx, this.prevVx, this.prevVx, this.prevVy, this.time);
   advect(2, this.Vy, this.prevVy, this.prevVx, this.prevVy, this.time);

   //clean up again
  project(this.Vx, this.Vy, this.prevVx, this.prevVy);

  diffuse(0, this.prevDensity, this.density, this.diffusion, this.time);
  advect(0,this.density, this.prevDensity, this.Vx, this.Vy, this.time); //ISSUEEEEEEEEE

  }
  
  void renderDensity() {
      for (int i = 0; i<N; i++) {
        for (int j = 0; j<N; j++) {
         int x = i * SCALE;
         int y = j * SCALE;
         float d = this.density[index(i,j)];         
         fill(d);
         noStroke();
         square(x,y,SCALE);  
        }
      }
  }
  
 void fadeDensity() {
   for (int i = 0; i< this.density.length; i++) {
     float d = density[i];
     density[i] = constrain(d-.1, 0, 255);
   }
 }
}
