#pragma once

#include <memory>
#include "Event.hpp"




namespace Engine {
	class Application {
	public:
		Application();
		virtual ~Application();

		Application(const Application&) = delete;
		Application(Application&&) = delete;
		Application& operator=(const Application&) = delete;
		Application& operator=(Application&&) = delete;


		virtual int start(unsigned int window_width, unsigned int window_height, const char* title);
		
		void close() { mbCloseWindow = true; };

		virtual void onUpdate() {
		};

		virtual void onUiDraw() {};

	private:
		void draw();
		std::shared_ptr<class Window> mpWindow;
		EventDispatcher mEventDispatcher;
		bool mbCloseWindow = false;

	};
}