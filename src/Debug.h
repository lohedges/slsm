// Zed Shaw's debug macros taken from "Learn C The Hard Way"
// http://c.learncodethehardway.org/book/ex20.html

#ifndef _DEBUG_H
#define _DEBUG_H

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>

/*! \file Debug.h
    \brief A set of debugging macros.
 */

#ifdef NDEBUG
#define slsm_debug(M, ...)
#else
#define slsm_debug(M, ...) fprintf(stderr, "[DEBUG] %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#endif

#define slsm_clean_errno() (errno == 0 ? "None" : strerror(errno))

#define slsm_log_err(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, slsm_clean_errno(), ##__VA_ARGS__)

#define slsm_log_warn(M, ...) fprintf(stderr, "[WARN] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, slsm_clean_errno(), ##__VA_ARGS__)

#define slsm_log_info(M, ...) fprintf(stderr, "[INFO] (%s:%d) " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)

#define slsm_check(A, M, ...) if(!(A)) { slsm_log_err(M, ##__VA_ARGS__); errno=0; goto error; }

#define slsm_sentinel(M, ...)  { slsm_log_err(M, ##__VA_ARGS__); errno=0; goto error; }

#define slsm_check_mem(A) check((A), "Out of memory.")

#define slsm_check_debug(A, M, ...) if(!(A)) { slsm_debug(M, ##__VA_ARGS__); errno=0; goto error; }

#endif /* _DEBUG_H */
