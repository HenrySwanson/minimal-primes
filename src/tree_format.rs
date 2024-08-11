pub struct Indenter {
    stack: Vec<Node>,
    prepared_pop: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Node {
    Live,
    Dead,
}

impl Indenter {
    const BAR: &'static str = " │  ";
    const TEE: &'static str = " ├─ ";
    const EMPTY: &'static str = "    ";
    const ELBOW: &'static str = " └─ ";

    pub fn new() -> Self {
        Self {
            stack: vec![],
            prepared_pop: false,
        }
    }

    pub fn write_line(&mut self, s: &str) {
        println!("{}{}", self.leader(), s);
        if self.prepared_pop {
            self.pop();
            self.prepared_pop = false;
        }
    }

    fn leader(&self) -> String {
        let mut s = String::new();
        for (idx, n) in self.stack.iter().copied().enumerate() {
            let is_last = idx == self.stack.len() - 1;
            s += match (n, is_last) {
                (Node::Live, true) => Self::TEE,
                (Node::Live, false) => Self::BAR,
                (Node::Dead, true) => Self::ELBOW,
                (Node::Dead, false) => Self::EMPTY,
            }
        }
        s
    }

    pub fn push_new(&mut self) {
        self.stack.push(Node::Live);
    }

    pub fn close(&mut self) {
        if let Some(last) = self.stack.last_mut() {
            *last = Node::Dead;
        }
    }

    pub fn pop(&mut self) {
        while let Some(last) = self.stack.last() {
            if *last == Node::Live {
                break;
            }
            self.stack.pop();
        }
    }

    pub fn prepare_pop(&mut self) {
        *self.stack.last_mut().expect("pop called when empty") = Node::Dead;
        self.close();
        assert!(!self.prepared_pop);
        self.prepared_pop = true;
    }
}
