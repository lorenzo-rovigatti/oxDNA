/*
 * PluginManager.h
 *
 *  Created on: 09/dic/2013
 *      Author: lorenzo
 */

#ifndef PLUGINMANAGER_H_
#define PLUGINMANAGER_H_

#include <string>
#include <vector>
#include <stack>

#include "../Observables/BaseObservable.h"
#include "../Utilities/parse_input/parse_input.h"
#include "../Interactions/BaseInteraction.h"
#include "../Backends/MCMoves/BaseMove.h"

/**
 * @brief Manages oxDNA plugins. As of now it only supports {@link BaseObservable observables} and {@link IBaseInteraction interactions}.
 * It implements the singleton pattern (see http://en.wikipedia.org/wiki/Singleton_pattern).
 *
 * As of now, only {@link BaseObservable observable}, {@link IBaseInteraction interaction} and
 * {@link BaseMove move}
 * plugins are supported. In order to write a plugin you should write a regular observable, interaction
 * or move class and add a simple function that serves as entry point for the plugin manager.
 * Note that the default name for the entry points is make_NAME, but make, make_observable, make_interaction
 * and make_move should also work (but be aware that they may cause name collisions).
 * As an example, we will assume that the new plugin is an observable named MyObservable. We write
 * this observable in two files, MyObservable.cpp and MyObservable.h. In order to provide the
 * required entry points we  add the following line at the end of the MyObservable.h file
 *
 @code
 extern "C" BaseObservable *make_MyObservable() { return new MyObservable(); }
 @endcode
 *
 * Plugins should be compiled as dynamic libraries. On linux systems and with the gcc
 * compiler this can be done as follows
 @code
 g++ -O3 -fpic -c -lm MyObservable.cpp -I/PATH/TO/OXDNA/src/Observables/
 g++ -shared -o MyObservable.so MyObservable.o
 @endcode
 * The final shared library, MyObservable.so, should be placed in a folder contained in the plugin
 * search path. This, by default, is the directory you run your simulation in. You can specify
 * additional paths by using the plugin_search_path key in the input file (see below for details).
 * In order to use the newly created observable in your simulation, use 'type = MyObservable' as
 * an observable specifier (as you would do for a regular observable). Note that the 'type' specifier
 * should match the case-sensitive filename and NOT the class name.
 *
 * The exact same procedure described above holds true for interactions and moves. For instance, the
 * entry function should be defined as follows (with the new class being named MyInteraction)
 *
 * @code
 extern "C" IBaseInteraction *make_MyInteraction() { return new MyInteraction(); }
 @endcode
 *
 * @verbatim
 [plugin_search_path = <string> (a semicolon-separated list of directories where plugins are looked for in, in addition to the current directory.)]
 [plugin_observable_entry_points = <string> (a semicolon-separated list of prefixes which will be used to look for entry points in shared libraries containing observables.)]
 [plugin_interaction_entry_points = <string> (a semicolon-separated list of prefixes which will be used to look for entry points in shared libraries containing interactions..)]
 [plugin_do_cleanup = <bool> (whether the PluginManager should perform a clean up at the end of the run. It should be set to false only when profiling code, since otherwise callgrind cannot read plugins' source files. It defaults to true.)]
 @endverbatim
 */
class PluginManager {
protected:
	static std::shared_ptr<PluginManager> _manager;
	std::vector<std::string> _path;
	bool _initialised;
	bool _do_cleanup;
	std::stack<void *> _handles;

	std::vector<std::string> _obs_entry_points;
	std::vector<std::string> _inter_entry_points;
	std::vector<std::string> _move_entry_points;

	void *_get_handle(std::string &name);
	void *_get_entry_point(void *handle, std::string name, std::vector<std::string> entry_points);

private:
	/**
	 * @brief The default constructor is private because this is a singleton class.
	 */
	PluginManager();
public:
	virtual ~PluginManager();

	PluginManager(PluginManager const &) = delete;

	/**
	 * @brief Initialises the singleton. Must be called before actually asking for plugins.
	 * @param sim_inp
	 */
	void init(input_file &sim_inp);

	/**
	 * @brief Add one or more directories to the plugin search path.
	 * @param s a semicolon-separated list of directories
	 */
	void add_to_path(std::string s);

	/**
	 * @brief Looks for an {@link BaseObservable observable} plugin in the current plugin path and, if found, builds it and returns it as a shared pointer
	 * @param name Case-sensitive name of the plugin
	 * @return a shared_ptr to the newly built plugin
	 */
	ObservablePtr get_observable(std::string name);

	/**
	 * @brief Looks for an {@link BaseInteraction interaction} plugin in the current plugin path and, if found, builds it and returns it as a shared pointer
	 * @param name Case-sensitive name of the plugin
	 * @return a shared_ptr to the newly built plugin
	 */
	InteractionPtr get_interaction(std::string name);

	/**
	 * @brief Looks for an {@link BaseInteraction interaction} plugin in the current plugin path and, if found, builds it and returns it as a shared pointer
	 * @param name Case-sensitive name of the plugin
	 * @return a shared_ptr to the newly built plugin
	 */
	MovePtr get_move(std::string name);

	/**
	 * @brief Returns a pointer to the current instance of PluginManager. Implements the singleton pattern.
	 * @return
	 */
	static std::shared_ptr<PluginManager> instance();
};

#endif /* PLUGINMANAGER_H_ */
