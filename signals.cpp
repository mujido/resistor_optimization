#include <atomic>
#include <csignal>

namespace
{
	std::atomic<bool> sigintFlag = false;

	extern "C" void sigintHandler(int)
	{
		sigintFlag.store(true, std::memory_order_release);
	}
}

void setupSignalHandler()
{
	signal(SIGINT, sigintHandler);
}

bool exitedRequested()
{
	return sigintFlag.load(std::memory_order_acquire);
}

