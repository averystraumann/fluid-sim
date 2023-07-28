
FluidCell fluid;


void settings() {
  size(N*SCALE, N*SCALE);
}



void setup() {
  fluid = new FluidCell(0.5, 0, 0.00001);
  
}

void mouseDragged() {
  fluid.addDye(mouseX/SCALE, mouseY/SCALE, 125);
  float amtX  = mouseX - pmouseX;
  float amtY  = mouseY - pmouseY;
  fluid.addVelocity(mouseX/SCALE, mouseY/SCALE,amtX,amtY);
  
}

void draw() {
  background(0);
  fluid.timeStep();
  fluid.renderDensity();
  fluid.fadeDensity();
  
}
