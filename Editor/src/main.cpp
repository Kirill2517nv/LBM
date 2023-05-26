#include <iostream>
#include <memory>
#include <Application.hpp>
#include <imgui/imgui.h>
#include <imgui_internal.h>
#include <implot/implot.h>
#include <vector>
#include <algorithm>

#pragma warning(disable: 6386)


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

		if (ImGui::BeginMenuBar())
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
		}
		ImGui::End();
	}

	void Draw_Heatmap() {

		std::vector<std::vector<std::vector<double>>> rho = Solver->get_rhomulticomponent();

		
		double scale_min = Solver->get_min_rho();
		double scale_max = Solver->get_max_rho();
		static int numspec = 0;

		static ImPlotColormap map = ImPlotColormap_Jet;
		if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), ImVec2(225, 0), map)) {
			map = (map + 1) % ImPlot::GetColormapCount();
		}
		ImGui::SameLine();
		ImGui::LabelText("##Colormap Index", "%s", "Change Colormap");
		
		

		static ImPlotAxisFlags axes_flags = ImPlotAxisFlags_Lock | ImPlotAxisFlags_NoGridLines | ImPlotAxisFlags_NoTickMarks;

		ImPlot::PushColormap(map);

		int size_x = Solver->get_Nx();
		int size_y = Solver->get_Ny();
		const int all_size = (size_x + 2) * (size_y + 2);
		double* values_rho = new double[all_size];
		for (int j = 0; j < size_y + 2; j++) {
			for (int i = 0; i < size_x + 2; i++) {
				values_rho[j * (size_x + 2) + i] = rho[numspec][i][j];
			}
		}

		if (ImPlot::BeginPlot("##Heatmap2", ImVec2(600, 300))) {
			ImPlot::SetupAxes(NULL, NULL, ImPlotAxisFlags_NoDecorations, ImPlotAxisFlags_NoDecorations);
			ImPlot::SetupAxesLimits(0, size_x + 1, 0, size_y + 1);
			ImPlot::PlotHeatmap("heat1", values_rho, size_y + 2, size_x + 2, 0, 0, NULL, ImPlotPoint(0, 0), ImPlotPoint(size_x + 1, size_y + 1));
			ImPlot::EndPlot();
		}
		ImGui::SameLine();
		ImPlot::ColormapScale("##HeatScale", scale_min, scale_max, ImVec2(60, 300));
		
		ImPlot::PopColormap();
		delete[] values_rho;

	}


	virtual void onUiDraw() override
	{
		double scale_min = Solver->get_min_rho();
		double scale_max = Solver->get_max_rho();
		setupDockspaceMenu();
		ImGui::Begin("My name is window, ImGUI window");
		ImGui::SetNextItemWidth(225);
		ImGui::InputDouble("min rho  ", &scale_min, 0.01f, 1.0f, "%.8f", ImGuiInputTextFlags_ReadOnly);
		ImGui::SameLine();
		ImGui::SetNextItemWidth(225);
		ImGui::InputDouble("max rho  ", &scale_max, 0.01f, 1.0f, "%.8f", ImGuiInputTextFlags_ReadOnly);
		ImGui::SameLine();
		ImGui::SetNextItemWidth(225);
		ImGui::InputInt("time", &time, 0.01f, 1.0f, ImGuiInputTextFlags_ReadOnly);
		Draw_Heatmap();
		ImGui::End();
	}
};

int main()
{
	srand(time(0));
	system("mkdir VTK");
	system("mkdir DATA");
	auto pEditor = std::make_shared<Editor>();
	int returnCode = pEditor->start(100, 5, 400, 2);

	return returnCode;
}