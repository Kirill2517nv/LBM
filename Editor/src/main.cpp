#include <iostream>
#include <memory>
#include <Application.hpp>
#include <imgui/imgui.h>
#include <imgui_internal.h>
#include <implot/implot.h>
#include <vector>


class Editor : public Engine::Application {

	virtual void onUpdate() override {

	}

	void setupDockspaceMenu()
	{
		static ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_PassthruCentralNode | ImGuiDockNodeFlags_NoWindowMenuButton;
		static ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;
		window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
		window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;
		window_flags |= ImGuiWindowFlags_NoBackground;

		const ImGuiViewport* viewport = ImGui::GetMainViewport();
		ImGui::SetNextWindowPos(viewport->WorkPos);
		ImGui::SetNextWindowSize(viewport->WorkSize);
		ImGui::SetNextWindowViewport(viewport->ID);
		ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
		ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
		ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
		ImGui::Begin("DockSpace", nullptr, window_flags);
		ImGui::PopStyleVar(3);

		ImGuiIO& io = ImGui::GetIO();
		ImGuiID dockspace_id = ImGui::GetID("MyDockSpace");
		ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);

		/*if (ImGui::BeginMenuBar())
		{
			if (ImGui::BeginMenu("File"))
			{
				if (ImGui::MenuItem("New Scene...", NULL))
				{

				}
				if (ImGui::MenuItem("Open Scene...", NULL))
				{

				}
				if (ImGui::MenuItem("Save Scene...", NULL))
				{

				}
				ImGui::Separator();
				if (ImGui::MenuItem("Exit", NULL))
				{
					close();
				}
				ImGui::EndMenu();
			}
			ImGui::EndMenuBar();
		}*/
		ImGui::End();
	}

	void Demo_Heatmaps() {
		std::vector<std::vector<double>> value = { {0.8f, 2.4f, 2.5f, 3.9f, 0.0f, 4.0f, 0.0f},
										{2.4f, 0.0f, 4.0f, 1.0f, 2.7f, 0.0f, 0.0f},
										{1.1f, 2.4f, 0.8f, 4.3f, 1.9f, 4.4f, 0.0f},
										{0.6f, 0.0f, 0.3f, 0.0f, 3.1f, 0.0f, 0.0f},
										{0.7f, 1.7f, 0.6f, 2.6f, 2.2f, 6.2f, 0.0f},
										{1.3f, 1.2f, 0.0f, 0.0f, 0.0f, 3.2f, 5.1f},
										{0.1f, 2.0f, 0.0f, 1.4f, 0.0f, 1.9f, 6.3f} };

		static float values1[7][7] = { {0.8f, 2.4f, 2.5f, 3.9f, 0.0f, 4.0f, 0.0f},
										{2.4f, 0.0f, 4.0f, 1.0f, 2.7f, 0.0f, 0.0f},
										{1.1f, 2.4f, 0.8f, 4.3f, 1.9f, 4.4f, 0.0f},
										{0.6f, 0.0f, 0.3f, 0.0f, 3.1f, 0.0f, 0.0f},
										{0.7f, 1.7f, 0.6f, 2.6f, 2.2f, 6.2f, 0.0f},
										{1.3f, 1.2f, 0.0f, 0.0f, 0.0f, 3.2f, 5.1f},
										{0.1f, 2.0f, 0.0f, 1.4f, 0.0f, 1.9f, 6.3f} };
		static float scale_min = 0;
		static float scale_max = 6.3f;
		static const char* xlabels[] = { "C1","C2","C3","C4","C5","C6","C7" };
		static const char* ylabels[] = { "R1","R2","R3","R4","R5","R6","R7" };

		static ImPlotColormap map = ImPlotColormap_Jet;
		if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), ImVec2(225, 0), map)) {
			map = (map + 1) % ImPlot::GetColormapCount();
		}

		ImGui::SameLine();
		ImGui::LabelText("##Colormap Index", "%s", "Change Colormap");
		ImGui::SetNextItemWidth(225);
		ImGui::DragFloatRange2("Min / Max", &scale_min, &scale_max, 0.01f, -20, 20);

		static ImPlotAxisFlags axes_flags = ImPlotAxisFlags_Lock | ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoTickMarks;

		ImPlot::PushColormap(map);

		const int size = 80;
		static double values2[size * size];
		srand((unsigned int)(ImGui::GetTime() * 1000000));
		for (int i = 0; i < size * size; ++i)
			values2[i] = rand() % 10 * 0.1;

		if (ImPlot::BeginPlot("##Heatmap2", ImVec2(600, 300))) {
			ImPlot::SetupAxes(NULL, NULL, ImPlotAxisFlags_NoDecorations, ImPlotAxisFlags_NoDecorations);
			ImPlot::SetupAxesLimits(0, 600, 0, 100);
			ImPlot::PlotHeatmap("heat1", values2, size, size, 0, 0, NULL, ImPlotPoint(0, 0), ImPlotPoint(600, 100));
			ImPlot::EndPlot();
		}
		ImGui::SameLine();
		ImPlot::ColormapScale("##HeatScale", scale_min, scale_max, ImVec2(60, 300));
		
		ImPlot::PopColormap();

	}


	virtual void onUiDraw() override
	{
		setupDockspaceMenu();
		ImGui::SetNextWindowSize(ImVec2(1600, 1600));
		ImGui::Begin("My name is window, ImGUI window");
		// Text that appears in the window
		ImGui::Text("Hello there adventurer!");
		Demo_Heatmaps();
		ImGui::End();
	}
};

int main()
{
	srand(time(0));
	system("mkdir VTK");
	system("mkdir DATA");
	auto pEditor = std::make_shared<Editor>();
	int returnCode = pEditor->start(1600, 900, "Editor");

	return returnCode;
}