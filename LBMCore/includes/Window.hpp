#pragma once
#include "Event.hpp"
#include <string>

#include <functional>

struct GLFWwindow;

namespace Engine {

	class Window {
	public:
		using EventCallbackFn = std::function<void(BaseEvent&)>;
		Window(std::string title, const unsigned int width, const unsigned int height);
		~Window();

		Window(const Window&) = delete;
		Window(Window&&) = delete;
		Window& operator=(const Window&) = delete;
		Window& operator=(Window&&) = delete;

		void onUpdate();
		unsigned int getWidth() const { return mData.width; };
		unsigned int getHeight() const { return mData.height; };

		void set_event_callback(const EventCallbackFn& callback) {
			mData.eventCallbackFn = callback;
		}


	private:
		struct WindowData {
			std::string title;
			unsigned int width;
			unsigned int height;
			EventCallbackFn eventCallbackFn;
		};
		int init();
		void shutdown();

		GLFWwindow* mpWindow = nullptr;
		WindowData mData;
	};
}