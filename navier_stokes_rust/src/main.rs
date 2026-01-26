mod grid;
mod renderer;

use grid::FluidGrid;
use renderer::FluidRenderer;
use piston_window::*;

fn main() {
    let mut grid = FluidGrid::new(40, 40, 15.0);
    // TODO: set some initial fluid motion.
    let center_i = grid.nx / 2;
    let center_j = grid.ny / 2;
    let idx = grid.v_idx(center_i, center_j);
    grid.v[idx] = 100.0;

    let mut viscosity: f64 = 0.0;
    let mut vorticity_strength: f64 = 5.0; // Default strength

    let mut window: PistonWindow = WindowSettings::new("Fluid Sim", [600, 600])
        .exit_on_esc(true)
        .resizable(true)
        .build()
        .unwrap();

    let mut show_speed = false;
    let mut mouse_down = false;
    let mut last_mouse_pos = [0.0, 0.0];

    while let Some(event) = window.next() {
        if let Some(Button::Keyboard(Key::C)) = event.press_args() {
            show_speed = !show_speed;
        }
        if let Some(Button::Keyboard(Key::R)) = event.press_args() {
            grid = FluidGrid::new(40, 40, 15.0); // Reset the grid to initial state.
            // Optionally, re-initialize any initial conditions here
            let center_i = grid.nx / 2;
            let center_j = grid.ny / 2;
            let idx = grid.v_idx(center_i, center_j);
            grid.v[idx] = 100.0;
        }
        if let Some(Button::Mouse(MouseButton::Left)) = event.press_args() {
            mouse_down = true;
        }
        if let Some(Button::Mouse(MouseButton::Left)) = event.release_args() {
            mouse_down = false;
        }
        if let Some(Button::Keyboard(Key::Right)) = event.press_args() {
            viscosity = (viscosity + 0.01).min(1.0_f64);
        }
        if let Some(Button::Keyboard(Key::Left)) = event.press_args() {
            viscosity = (viscosity - 0.01).max(0.0_f64);
        }
        if let Some(Button::Keyboard(Key::Up)) = event.press_args() {
            vorticity_strength += 1.0;
        }
        if let Some(Button::Keyboard(Key::Down)) = event.press_args() {
            vorticity_strength = (vorticity_strength - 1.0).max(0.0);
        }
        if let Some([x, y]) = event.mouse_cursor_args() {
            if mouse_down {
                // Inject velocity based on mouse movement
                let dx = x - last_mouse_pos[0];
                let dy = y - last_mouse_pos[1];
                grid.add_velocity_at_pixel(x, y, dx, dy);
            }
            last_mouse_pos = [x, y];
        }
        if let Some(_args) = event.update_args() {
            grid.run_step(0.016, viscosity, vorticity_strength);
        }

        window.draw_2d(&event, |context, graphics, _device| {
            clear([0.1, 0.1, 0.1, 1.0], graphics); // Clear to dark gray

            // We'll create and call our drawing function here
            if show_speed {
                FluidRenderer::draw_speed(&grid, &context, graphics);
            } else {
                FluidRenderer::draw_velocities(&grid, &context, graphics);
            }
        });
    }

}
