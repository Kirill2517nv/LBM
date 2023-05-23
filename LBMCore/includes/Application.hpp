#pragma once

#include <memory>
#include "Event.hpp"
#include "BasicSolver2D.hpp"



namespace Engine {
	class Application {
	public:
		Application();
		virtual ~Application();

		Application(const Application&) = delete;
		Application(Application&&) = delete;
		Application& operator=(const Application&) = delete;
		Application& operator=(Application&&) = delete;


		virtual int start(int Nx, int Ny, double T, int numspec, 
			unsigned int window_width = 1600, unsigned int window_height = 900, const char* title = "Solver");
		
		void close() { mbCloseWindow = true; };

		virtual void onUpdate() {
		};

		virtual void onUiDraw() {};

	private:
		void draw();
		std::shared_ptr<class Window> mpWindow;
		EventDispatcher mEventDispatcher;
		bool mbCloseWindow = false;
	protected:
		std::shared_ptr<BasicSolver2D> Solver;
		double* rho;
	};
}