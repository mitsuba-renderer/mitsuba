static inline float generateRandomNumber()
{
	const float U = ((float)rand()) / (float)RAND_MAX;
	return U;
}

static inline bool IsFiniteNumber(float x) 
{
	return (x <= std::numeric_limits<Float>::max() && x >= -std::numeric_limits<Float>::max()); 
} 


static inline double  abgam (double x)
{
  double  gam[10],
          temp;

  gam[0] = 1./ 12.;
  gam[1] = 1./ 30.;
  gam[2] = 53./ 210.;
  gam[3] = 195./ 371.;
  gam[4] = 22999./ 22737.;
  gam[5] = 29944523./ 19733142.;
  gam[6] = 109535241009./ 48264275462.;
  temp = 0.5*log (2*M_PI) - x + (x - 0.5)*log (x)
    + gam[0]/(x + gam[1]/(x + gam[2]/(x + gam[3]/(x + gam[4] /
	  (x + gam[5]/(x + gam[6]/x))))));

  return temp;
}

static inline double  gamma (double x)
{
  double  result;
  result = exp (abgam (x + 5))/(x*(x + 1)*(x + 2)*(x + 3)*(x + 4));
  return result;
}

static inline double  beta (double m, double n)
{
  return (gamma (m)*gamma (n)/gamma (m + n));
}

#define vec3 Vector
#define vec2 Vector

struct RayInfo
{
	// direction
	vec3 w;
	float theta;
	float cosTheta;
	float sinTheta;
	float tanTheta;
	float alpha;
	float Lambda;

	void updateDirection(const vec3& w, const float alpha_x, const float alpha_y)
	{
		this->w = w;
		theta = acosf(w.z);
		cosTheta = w.z;
		sinTheta = sinf(theta);
		tanTheta = sinTheta / cosTheta;
		const float invSinTheta2 = 1.0f / (1.0f - w.z*w.z);
		const float cosPhi2 = w.x*w.x*invSinTheta2;
		const float sinPhi2 = w.y*w.y*invSinTheta2;
		alpha = sqrtf( cosPhi2*alpha_x*alpha_x + sinPhi2*alpha_y*alpha_y ); 
		// Lambda
		if(w.z > 0.9999f)
			Lambda = 0.0f;
		else if(w.z < -0.9999f)
			Lambda = -1.0f;
		else
		{	
			const float a = 1.0f/tanTheta/alpha;
			Lambda = 0.5f*(-1.0f + ((a>0)?1.0f:-1.0f) * sqrtf(1 + 1/(a*a)));
		}
	}

	// height
	float h;
	float C1;
	float G1;

	void updateHeight(const float& h)
	{
		this->h = h;
		C1 = std::min(1.0f, std::max(0.0f, 0.5f*(h+1.0f)));

		if(this->w.z > 0.9999f)
			G1 = 1.0f;
		else if(this->w.z <= 0.0f)
			G1 = 0.0f;
		else
			G1 = powf(this->C1, this->Lambda);
	}
};

inline float invC1(const float U) 
{
	const float h = std::max(-1.0f, std::min(1.0f, 2.0f*U-1.0f));
	return h;	
}

inline float sampleHeight(const RayInfo& ray, const float U) 
{
	if(ray.w.z > 0.9999f)
		return std::numeric_limits<Float>::max();
	if(ray.w.z < -0.9999f)
	{
		const float value = invC1(U*ray.C1);
		return value;
	}
	if(fabsf(ray.w.z) < 0.0001f)
		return ray.h;

	// probability of intersection
	if (U > 1.0f - ray.G1) // leave the microsurface
		return std::numeric_limits<Float>::max();

	const float h = invC1( 
			ray.C1 / powf((1.0f-U),1.0f/ray.Lambda)
			);
	return h;
}

float D_ggx(const vec3& wm, const float alpha_x, const float alpha_y) {
	if( wm.z <= 0.0f)
		return 0.0f;

	// slope of wm
	const float slope_x = -wm.x/wm.z;
	const float slope_y = -wm.y/wm.z;

	// P22
	const float tmp = 1.0f + slope_x*slope_x/(alpha_x*alpha_x) + slope_y*slope_y/(alpha_y*alpha_y);
	const float P22 = 1.0f / (M_PI * alpha_x * alpha_y) / (tmp * tmp);

	// value
	const float value = P22 / (wm.z*wm.z*wm.z*wm.z);
	return value;
}

vec2 sampleP22_11(const float theta_i, const float U, const float U_2, const float alpha_x, const float alpha_y)
{
	vec2 slope;

	if(theta_i < 0.0001f)
	{
		const float r = sqrtf(U/(1.0f-U));
		const float phi = 6.28318530718f * U_2;
		slope.x = r * cosf(phi);
		slope.y = r * sinf(phi);
		return slope;
	}

	// constant
	const float sin_theta_i = sinf(theta_i);
	const float cos_theta_i = cosf(theta_i);
	const float tan_theta_i = sin_theta_i/cos_theta_i;

	// slope associated to theta_i
	const float slope_i = cos_theta_i/sin_theta_i;

	// projected area
	const float projectedarea = 0.5f * (cos_theta_i + 1.0f);
	if(projectedarea < 0.0001f || projectedarea!=projectedarea)
		return vec2(0,0,0);
	// normalization coefficient
	const float c = 1.0f / projectedarea;

	const float A = 2.0f*U/cos_theta_i/c - 1.0f;
	const float B = tan_theta_i;
	const float tmp = 1.0f / (A*A-1.0f);

	const float D = sqrtf(std::max(0.0f, B*B*tmp*tmp - (A*A-B*B)*tmp));
	const float slope_x_1 = B*tmp - D;
	const float slope_x_2 = B*tmp + D;
	slope.x = (A < 0.0f || slope_x_2 > 1.0f/tan_theta_i) ? slope_x_1 : slope_x_2;

	float U2;
	float S;
	if(U_2 > 0.5f)
	{
	S = 1.0f;
	U2 = 2.0f*(U_2-0.5f);
	}
	else
	{
	S = -1.0f;
	U2 = 2.0f*(0.5f-U_2);
	}
	const float z = (U2*(U2*(U2*0.27385f-0.73369f)+0.46341f)) / (U2*(U2*(U2*0.093073f+0.309420f)-1.000000f)+0.597999f);
	slope.y = S * z * sqrtf(1.0f+slope.x*slope.x);

	return slope;
}

vec3 sampleVNDF(const vec3& wi, const float alpha_x, const float alpha_y)
{
	const float U1 = generateRandomNumber();
	const float U2 = generateRandomNumber();

	// sample D_wi

	// stretch to match configuration with alpha=1.0	
	const vec3 wi_11 = normalize(vec3(alpha_x * wi.x, alpha_y * wi.y, wi.z));

	// sample visible slope with alpha=1.0
	vec2 slope_11 = sampleP22_11(acosf(wi_11.z), U1, U2, alpha_x, alpha_y);

	// align with view direction
	const float phi = atan2(wi_11.y, wi_11.x);
	vec2 slope(cosf(phi)*slope_11.x - sinf(phi)*slope_11.y, sinf(phi)*slope_11.x + cos(phi)*slope_11.y, 0);

	// stretch back
	slope.x *= alpha_x;
	slope.y *= alpha_y;

	// if numerical instability
	if( (slope.x != slope.x) || !IsFiniteNumber(slope.x) ) 
	{
		if(wi.z > 0) return vec3(0.0f,0.0f,1.0f);
		else return normalize(vec3(wi.x, wi.y, 0.0f));
	}

	// compute normal
	const vec3 wm = normalize(vec3(-slope.x, -slope.y, 1.0f));

	return wm;
}
