/**
 * @file 	SDLDevice.h
 *
 * @date 	23.2.2013
 * @author 	Zdenek Travnicek <travnicek@iim.cz>
 * @copyright GNU Public License 3.0
 *
 */

#ifndef SDLDEVICE_H_
#define SDLDEVICE_H_
#include "PlatformDefs.h"
#include "keys.h"
#include <string>
#include <thread>
#include <mutex>
#include <atomic>
#include <vector>
#include <memory>

namespace iimavlib {

struct EXPORT RGB {
	uint8_t r:8;
	uint8_t g:8;
	uint8_t b:8;
};
struct sdl_pimpl_t;
class EXPORT SDLDevice {
public:
	typedef std::vector<RGB> data_type;
	SDLDevice(size_t width, size_t height, const std::string& title = "IIMAudio application", bool fullscreen = false);
	virtual ~SDLDevice();
	bool start();
	bool stop();
	template<typename T>
	bool update(const std::vector<T>&);
	bool is_stopped() const;
private:
	size_t width_;
	size_t height_;
	const std::string title_;
	std::thread thread_;
	std::mutex thread_mutex_;
	std::mutex surface_mutex_;
	std::mutex data_mutex_;
	std::atomic<bool> finish_;
	std::unique_ptr<sdl_pimpl_t> pimpl_;
	data_type data_;
	bool data_changed_;
	bool fullscreen_;
	std::atomic<bool> flip_required_;

	void run();
	bool process_events();
	void update_data();
	bool key_pressed(const int key, bool pressed);
	bool mouse_moved(const int x, const int y, const int dx, const int dy);
	bool mouse_button(const int key, const bool pressed, const int x, const int y);
	virtual bool do_key_pressed(const int key, bool pressed);
	virtual bool do_mouse_moved(const int x, const int y, const int dx, const int dy);
	virtual bool do_mouse_button(const int key, const bool pressed, const int x, const int y);
};

template<>
EXPORT bool SDLDevice::update(const data_type& data);

template<typename T>
bool SDLDevice::update(const std::vector<T>& data) {
	if (finish_) return false;
	std::unique_lock<std::mutex> lock(data_mutex_);
	std::copy(reinterpret_cast<data_type::pointer>(&data[0]),
			reinterpret_cast<data_type::pointer>(&data[0])+data.size()*sizeof(T)/sizeof(RGB),
			&data_[0]);
	data_changed_ = true;
	return true;
}

}


#endif /* SDLDEVICE_H_ */

