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

#include "../Observables/BaseObservable.h"
#include "../Utilities/parse_input/parse_input.h"

/**
 * @brief Manages oxDNA plugins. As of now it only supports {@link BaseObservable observables}.
 * It implements the singleton pattern (see http://en.wikipedia.org/wiki/Singleton_pattern).
 *
 * As of now, only {@link BaseObservable observables} plugins are supported. In order to write an
 * observable plugin you should write a regular observable and add two simple functions that serve as
 * entry points for the plugin manager.
 * As an example, we will assume that the new observable plugin is called MyObservable. We write
 * this observable in two files, MyObservable.cpp and MyObservable.h. In order to provide the
 * required entry points we  add the following two lines at the end of the MyObservable.h file
 *
@code
extern "C" BaseObservable<float> *make_float() { return new MyObservable<float>(); }
extern "C" BaseObservable<double> *make_double() { return new MyObservable<double>(); }
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
 * @verbatim
[plugin_search_path = <string> (a semicolon-separated list of directories where plugins are looked for in, in addition to the current directory.)]
@endverbatim
 */
class PluginManager {
protected:
	static PluginManager *_manager;
	std::vector<std::string> _path;
	bool _initialised;

private:
	PluginManager();
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
