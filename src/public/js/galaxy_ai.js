import * as THREE from 'https://unpkg.com/three/build/three.module.js';

class Star {
    constructor(x, y, z, type) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.type = type;
    }

    getColor() {
        const colors = {
            'O': 0xFF0000, // Красный
            'B': 0x00FFFF, // Синий
            'A': 0x00FF00, // Зелёный
            'F': 0xFFFF00, // Жёлтый
            'G': 0xFFFFFF, // Белый
            'K': 0xFFA500, // Оранжевый
            'M': 0x800000 // Тёмно-красный
        };
        return colors[this.type] || 0xFFFFFF;
    }
}

class GalaxyGenerator {
    constructor(seed = 12345) {
        this.seed = seed;
    }

    // Linear Congruential Generator (LCG) for pseudorandom numbers
    lcg(seed) {
        return () => {
            seed = (214013 * seed + 2531011) % 2147483648;
            return seed / 2147483648;
        };
    }

    // Beta distribution generator (rejection method)
    betaDistribution(alpha, beta, randFunc) {
        let x, y;
        do {
            x = Math.pow(randFunc(), 1 / alpha);
            y = Math.pow(randFunc(), 1 / beta);
        } while (x + y > 1.0);

        return x / (x + y);
    }

    hashCode(str) {
        let hash = 0;
        for (let i = 0; i < str.length; i++) {
            hash = ((hash << 5) - hash) + str.charCodeAt(i);
            hash |= 0; // Convert to 32bit integer
        }
        return hash;
    }

    generateStarPatch(x0, y0, z0, t0, distance = 100) {
        const sectorSize = 100;
        const sectorX = Math.floor(x0 / sectorSize);
        const sectorY = Math.floor(y0 / sectorSize);
        const sectorZ = Math.floor(z0 / sectorSize);

        // Pseudo-random generator initialization
        const seed = this.hashCode([sectorX, sectorY, sectorZ].join(',')) % 2147483648;
        const rand = this.lcg(seed);

        // Parameters generation
        const N = Math.floor(-50 * Math.log(1 - rand()));
        const stars = [];

        const G = 1; // Gravitational constant in normalized units
        const Omega0 = 0.2;
        const r0 = 1500;

        for (let i = 0; i < N; i++) {
            // Orbital parameters
            const a = -1500 * Math.log(1 - rand());
            const e = this.betaDistribution(2, 5, rand);
            const inc = Math.acos(2 * rand() - 1);
            const omega = 2 * Math.PI * rand();

            // Current position calculation
            const r = a * (1 - e ** 2) / (1 + e * Math.cos(omega));
            const Omega_p = Omega0 * Math.exp(-r / r0);
            const theta = omega + Omega_p * t0 + Math.sqrt(G * 1 / a ** 3) * t0;

            // XYZ coordinates
            const x = r * Math.cos(theta) + x0;
            const y = r * Math.sin(theta) + y0;
            const z = a * Math.sin(inc) * Math.sin(omega + Omega_p * t0) + z0;

            // Distance check
            if (Math.sqrt((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2) <= distance) {
                const type = this.getStarType(r);
                stars.push(new Star(x, y, z, type));
            }
        }

        return stars;
    }

    getStarType(r) {
        // Simple classification by distance from galaxy center.
        const types = ['O', 'B', 'A', 'F', 'G', 'K', 'M'];
        const index = Math.min(Math.floor(r / 500), types.length - 1);
        return types[index];
    }
}

// Visualization with Three.js
function initThreeJS(stars) {
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(30, window.innerWidth / window.innerHeight, 0.1, 1000);
    const renderer = new THREE.WebGLRenderer();

    renderer.setSize(window.innerWidth, window.innerHeight);
    document.body.appendChild(renderer.domElement);

    stars.forEach(star => {
        const geometry = new THREE.SphereGeometry(.1, 8, 8);
        const material = new THREE.MeshBasicMaterial({color: star.getColor()});
        const sphere = new THREE.Mesh(geometry, material);
        sphere.position.set(star.x / 100, star.y / 100, star.z / 100);
        scene.add(sphere);
    });

    camera.position.x = 0;
    camera.position.z = 200;

//    let rX = 0, drX = Math.PI / 1800;
    function animate() {
//	camera.rotation.x = rX;
//	camera.updateProjectionMatrix();
//	rX -= drX;
//	if (rX < 0) {
//	    rX = Math.PI * 2;
//	}
        requestAnimationFrame(animate);
        renderer.render(scene, camera);
    }

    animate();
}

/**
 * @typedef {Object} XY
 * @property {number} x
 * @property {number} y
 */


/**
 * @param {number} n
 * @return XY[]
 */
function generateAngles(n)
{
    let angle = Math.PI * 2 / n;
    const angles = [];
    for (let i = 0; i < n; ++i) {
        angles.push({x: Math.cos(angle * i), y: Math.sin(angle * i)});
    }

    return angles;
}

const generator = new GalaxyGenerator();
const time = 1e9;
const R = 400;
let stars = generator.generateStarPatch(0, 0, 0, time, R);
generateAngles(6).forEach(xy => stars = stars.concat(generator.generateStarPatch(R * xy.x, R * xy.y, 0, time, R)));
console.log(stars);
initThreeJS(stars);
