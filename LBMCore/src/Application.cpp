#include <iostream>
#include <fstream>


#include "Application.hpp"
#include "Log.hpp"
#include "Window.hpp"
#include "UIModule.hpp"
#include "Event.hpp"


#include <glad/glad.h>
#include <imgui/imgui.h>
#include <GLFW/glfw3.h>



namespace Engine {
	
	const float textScaleS = 10;

	Application::Application() {
		LOG_INFO("Starting application");
	}

	Application::~Application() {
		LOG_INFO("Closing application");
	}

	int Application::start(unsigned int window_width, unsigned int window_height, const char* title)
	{	
		// making a window
		mpWindow = std::make_shared<Window>(title, window_width, window_height);


		// Events ( add new event callback here and event itself in Event class)
		mEventDispatcher.addEventListener<EventMouseMoved>(
			[](EventMouseMoved& event) {
				//LOG_INFO("[MouseMoved] Mouse moved to {0} x {1}", event.x, event.y);
			}
		);

		mEventDispatcher.addEventListener<EventWindowResize>(
			[&](EventWindowResize& event) {
				LOG_INFO("[WindowResized] Changed size to {0} x {1}", event.width, event.height);
				draw();
			}
		);

		mEventDispatcher.addEventListener<EventWindowClose>(
			[&](EventWindowClose& event) {
				LOG_INFO("[WindowClose]");
				close();
			}
		);

		mpWindow->set_event_callback(
			[&](BaseEvent& event) {
				mEventDispatcher.dispatch(event);
			}
		);
		
		// main cycle
		while (!mbCloseWindow) {
			draw();
		}
		mpWindow = nullptr;

		return 0;
	}

	
	void Application::draw() {

		glClearColor(0.33f, 0.33f, 0.33f, 0);
		glClear(GL_COLOR_BUFFER_BIT);
		UIModule::onUiDrawBegin();
		onUiDraw();
		UIModule::onUiDrawEnd();
		mpWindow->onUpdate();
		onUpdate();
		

	}
}