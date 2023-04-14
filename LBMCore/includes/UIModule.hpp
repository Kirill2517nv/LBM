#pragma once

struct GLFWwindow;

namespace Engine {

    class UIModule {
    public:
        static void onWindowCreate(GLFWwindow* pWindow);
        static void onWindowClose();
        static void onUiDrawBegin();
        static void onUiDrawEnd();
    };

}