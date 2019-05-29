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

/**
 * @brief Manages oxDNA plugins. As of now it only supports {@link BaseObservable observables} and {@link IBaseInteraction interactions}.
 * It implements the singleton pattern (see http://en.wikipedia.org/wiki/Singleton_pattern).
 *
 * As of now, only {@link BaseObservable observable} and {@link IBaseInteraction interaction}
 * plugins are supported. In order to write a plugin you should write a regular observable
 * or interaction class and add two simple functions that serve as entry points for the plugin
 * manager. Note that the default names for the entry points are make_NAME_\<precision\> for
 * both observables and interactions or make_\<precision\> and make_observable_\<precision\> for
 * observables and make_\<precision\> and make_interaction\<precision\> for the interactions,
 * where \<precision\> should be either float or double.
 * As an example, we will assume that the new plugin is an observable named MyObservable. We write
 * this observable in two files, MyObservable.cpp and MyObservable.h. In order to provide the
 * required entry points we  add the following two lines at the end of the MyObservable.h file
 *
@code
extern "C" BaseObservable<float> *make_MyObservable_float() { return new MyObservable<float>(); }
extern "C" BaseObservable<double> *make_MyObservable_double() { return new MyObservable<double>(); }
@endcode
 *
 * Observable plugins should be compiled as dynamic libraries. On linux systems and with the gcc
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
 * The exact same procedure described above holds true for interactions. In this case, the two
 * entry functions should be defined as follows (with the new class being named MyInteraction)
 *
 * @code
extern "C" IBaseInteraction<float> *make_float() { return new MyInteraction<float>(); }
extern "C" IBaseInteraction<double> *make_double() { return new MyInteraction<double>(); }
@endcode
 *
 * @verbatim
[plugin_search_path = <string> (a semicolon-separated list of directories where plugins are looked for in, in addition to the current directory.)]
[plugin_observable_entry_points = <string> (a semicolon-separated list of prefixes which will be used to look for entry points in shared libraries containing observables, followed by either float or double.)]
[plugin_interaction_entry_points = <string> (a semicolon-separated list of prefixes which will be used to look for entry points in shared libraries containing interactions, followed by either float or double.)]
[plugin_do_cleanup = <bool> (whether the PluginManager should perform a clean up at the end of the run. It should be set to false only when profiling code, since otherwise callgrind cannot read plugins' source files. It defaults to true.)]
@endverbatim
 */
class PluginManager {
protected:
	static PluginManager *_manager;
	std::vector<std::string> _path;
	bool _initialised;
	bool _do_cleanup;
	std::stack<void *> _handles;

	std::vector<std::string> _obs_entry_points;
	std::vector<std::string> _inter_entry_points;

	void *_get_handle(std::string &name);
	void *_get_entry_point(void *handle, std::string name, std::vector<std::string> entry_points, std::string suffix);

private:
	/**
	 * @brief The default constructor is private because this is a singleton class.
	 */
	PluginManager();
	/**
	 * @brief The copy constructor is private because this is a singleton class.
	 * @param
	 */
	PluginManager(PluginManager const &);
public:
	virtual ~PluginManager();

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
	 * @brief Looks for an {@link BaseObservable observable} plugin in the current plugin path and, if found, builds it and returns it as a pointer
	 * @param name Case-sensitive name of the plugin
	 * @return a pointer to the newly built plugin
	 */
	template<typename number>
	BaseObservable<number> *get_observable(std::string name);

	/**
	 * @brief Looks for an {@link BaseInteraction interaction} plugin in the current plugin path and, if found, builds it and returns it as a pointer
	 * @param name Case-sensitive name of the plugin
	 * @return a pointer to the newly built plugin
	 */
	template<typename number>
	IBaseInteraction<number> *get_interaction(std::string name);

	/**
	 * @brief Cleans up the manager.
	 */
	static void clear();

	/**
	 * @brief Returns a pointer to the current instance of PluginManager. Implements the singleton pattern.
	 * @return
	 */
	static PluginManager *instance();
};

#endif /* PLUGINMANAGER_H_ */
